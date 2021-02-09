//
// Created by pablo on 30/10/20.
//

#include "telomereMutations.h"

#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)




int
cmri::mainTelomereMutations(common_options_t common_options, telomere_mutations_options_t telomere_mutation_options) {

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 3;

    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment
    mopt.flag |= MM_F_OUT_CG; // perform alignment


    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile f = gzopen(telomere_mutation_options.query_file.c_str(), "r");
    assert(f);
    kseq_t *ks = kseq_init(f);

    std::vector<mutations_t> mutations;

    // open index reader
    mm_idx_reader_t *index_reader = mm_idx_reader_open(telomere_mutation_options.target_file.c_str(), &iopt, 0);
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(index_reader, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt,
                         mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
        gzrewind(f);
        kseq_rewind(ks);
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            mm_reg1_t *reg;
            int j, i, n_reg;
            reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query

            mutations_t mut;
            mut.score=0;

            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                mm_reg1_t *r = &reg[j];
                assert(r->p); // with MM_F_CIGAR, this should not be NULL

                LOGGER.debug << ks->name.s << "\tlen:"
                            << ks->seq.l << "\tqs:"
                            << r->qs << "\tqe:"
                            << r->qe << "\t"
                            << "+-"[r->rev] << std::endl;

                LOGGER.debug << mi->seq[r->rid].name << "\tlen:"
                            << mi->seq[r->rid].len << "\trs:"
                            << r->rs << "\tre:"
                            << r->re << "\tmlen:"
                            << r->mlen << "\tblen:"
                            << r->blen << "\tmapq:"
                            << r->mapq << "\t cg:";
                for (i = 0; i <
                            r->p->n_cigar; ++i) {// IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                    LOGGER.debug << (r->p->cigar[i] >> 4) << "MIDNSH"[r->p->cigar[i] & 0xf];
                }
                LOGGER.debug << std::endl;

                void *km = nullptr;
                char *cs_str = NULL;
                int max_len = 0;
                int n_cs = mm_gen_cs(km, &cs_str, &max_len, mi, r, ks->seq.s, true);

                LOGGER.debug << ks->seq.s << std::endl;
                LOGGER.debug << ks->qual.s << std::endl;
                LOGGER.debug << cs_str << std::endl;


                auto cm_tag = parseTag(cs_str);
                std::map<int, std::vector<sbs_t>> sbs;
                std::map<int, std::vector<indel_t>> indels;

                int rs = r->rs;
                int re = r->re;
                int qs = r->qs;
                int qe = r->qe;
                std::string que;
                std::vector<int> qv;
                if(r->rev){
                    qs = ks->seq.l - r->qe;
                    qe = ks->seq.l - r->qs;
                    for (int i = ks->qual.l-1; i >= 0 ; i--) {
                        qv.push_back(static_cast<int>(ks->qual.s[i]) - 33);
                    }
                    que=reverse_complement(ks->seq.s);
                }
                else{
                    for (int i = 0; i < ks->qual.l; i++) {
                        qv.push_back(static_cast<int>(ks->qual.s[i]) - 33);
                    }
                    que = ks->seq.s;
                }
                int size;
                for (auto &item : cm_tag) {
                    switch (item.first) {
                        case ':':
                            size = std::stoi(item.second);
                            LOGGER.debug << "MATCH: " << que.substr(qs, size) << std::endl;
                            rs += size;
                            qs += size;
                            break;
                        case '-':
                        case '+': {
                            size = item.second.size();
                            double sum = 0;
                            for (int i = qs; i < qs + size; i++) { sum += qv[i]; }
                            std::string kind = item.first == '-' ? "del" : "ins";
                            LOGGER.debug << kind << " " << item.second << std::endl;
                            indel_t del = {kind, rs % 6, item.second, sum / size};
                            indels[int(std::ceil(rs / 6))].push_back(del);
                            if (item.first == '-') {
                                rs += size;
                            } else {
                                qs += size;
                            }
                        }
                            break;
                        case '*': {
                            LOGGER.debug << "SBS: " << item.second << " " << que[qs] << std::endl;
                            sbs_t m = {rs % 6, que[qs], qv[qs]};
                            sbs[int(std::ceil(rs / 6))].push_back(m);
                            size = item.second.size() / 2;
                            rs += size;
                            qs += size;
                        }
                            break;
                        default:
                            LOGGER.error << "Unknown state: " << item.first;
                            exit(-1);
                    }

                }

                double score = ( static_cast<double>(r->blen+ r->mlen)/ks->seq.l + r->mapq/60.0)/3.0;
                if(mut.score < score) {
                    mut.name = ks->name.s;
                    mut.rs = r->rs;
                    mut.re = r->re;
                    mut.qs = r->qs;
                    mut.qe = r->qe;
                    mut.seq_len = ks->seq.l;
                    mut.mapq = r->mapq;
                    mut.indels = indels;
                    mut.sbs = sbs;
                    mut.blen = r->blen;
                    mut.mlen = r->mlen;
                    mut.score = score;
                    mut.variants = find_variants(que);
                    mut.seq = que;
                    mut.cs_str = cs_str;
                }

                free(cs_str);
                free(r->p);
            }

            if(mut.score > 0) {
                mut.find_mutations();
                mutations.push_back(mut);
            }

            free(reg);
        }


        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(mi);
    }


    std::ofstream file(common_options.output_path + "/output.json");
    file << cmri::serialize(mutations);
    file.close();

    mm_idx_reader_close(index_reader); // close the index reader
    kseq_destroy(ks); // close the query file
    gzclose(f);
    return 0;


}
