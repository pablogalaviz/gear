//
// Created by pablo on 30/10/20.
//

#include "telomereMutations.h"

#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
#include <cmath>

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
    std::vector<char> wt_motif;
    for(auto c : telomere_mutation_options.wt_motif){
        wt_motif.push_back(c);
    }
    int wt_size = wt_motif.size();

    // open index reader
    mm_idx_reader_t *index_reader = mm_idx_reader_open(telomere_mutation_options.target_file.c_str(), &iopt, 0);
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(index_reader, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt,
                         mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
        gzrewind(f);
        kseq_rewind(ks);
        int count=0;
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            count++;
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


                auto cm_tag = parseTag(cs_str,r->rev);
                std::map<int, std::vector<sbs_t>> sbs;
                std::map<int, std::vector<indel_t>> indels;

                int rs = r->rs;
                int re = r->re;
                int qs = r->qs;
                int qe = r->qe;
                std::string que;
                std::string qv_str;
                std::vector<int> qv;
                if(r->rev){
                    qs = ks->seq.l - r->qe;
                    qe = ks->seq.l - r->qs;
                    for (int i = ks->qual.l-1; i >= 0 ; i--) {
                        qv.push_back(static_cast<int>(ks->qual.s[i]) - 33);
                    }
                    que=reverse_complement(ks->seq.s);
                    qv_str = reverse_string(ks->qual.s);
                }
                else{
                    for (int i = 0; i < ks->qual.l; i++) {
                        qv.push_back(static_cast<int>(ks->qual.s[i]) - 33);
                    }
                    que = ks->seq.s;
                    qv_str = ks->qual.s;
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
                            indel_t indel;
                            indel.kind = kind;
                            indel.pos = rs%wt_size;
                            indel.seq = item.second;
                            indel.mean_qv = sum/size;
                            indels[qs].push_back(indel);
                            if (item.first == '-') {
                                rs += size;
                            } else {
                                qs += size;
                            }
                        }
                            break;
                        case '*': {
                            LOGGER.debug << "SBS: " << item.second << " " << que[qs] << std::endl;
                            int is=std::max<int>(0,qs-2);
                            int ie=std::min<int>(qv.size(),qs+3);
                            double sum=0;
                            int total=0;
                            for (int i = is; i < ie; i++) { sum += qv[i];total++;}
                            sbs_t m;
                            m.pos=rs % wt_size;
                            m.value = que[qs];
                            m.qv = qv[qs];
                            m.mean_qv = sum / total;
                            sbs[qs].push_back(m);
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

                    auto rng = get_trimmed_range(qv,telomere_mutation_options.trimming_window_mean,telomere_mutation_options.trimming_threshold);

                    for (auto &item : sbs)
                        for (auto &v : item.second) {
                            if (item.first < rng.first || item.first > rng.second) { v.is_trimmed = true; }
                        }
                    for (auto &item : indels)
                        for (auto &v : item.second) {
                            if (item.first < rng.first || item.first > rng.second) { v.is_trimmed = true; }
                        }

                    mut.name = ks->name.l>0?ks->name.s:"";
                    mut.comment = ks->comment.l >0 ? ks->comment.s : "";
                    mut.rs = r->rs;
                    mut.re = r->re;
                    mut.qs = r->rev ? ks->seq.l - r->qe : r->qs;
                    mut.qe = r->rev ? ks->seq.l - r->qs : r->qe;
                    mut.ts = rng.first;
                    mut.te = rng.second;
                    mut.seq_len = ks->seq.l;
                    mut.mapq = r->mapq;
                    mut.indels = indels;
                    mut.sbs = sbs;
                    mut.blen = r->blen;
                    mut.mlen = r->mlen;
                    mut.score = score;
                    mut.variants = find_variants(que.substr(mut.qs,mut.qe-mut.qs));
                    mut.seq = que;
                    mut.qv = qv_str;
                    mut.seq_trimmed=trimm_string(que,rng.first,rng.second,'N');
                    mut.qv_trimmed=trimm_string(qv_str,rng.first,rng.second,'!');
                    mut.cs_str = cs_str;
                    mut.reverse = r->rev;
                }

                free(cs_str);
                free(r->p);
            }

            if(mut.score > 0) {
                mut.find_mutations(wt_motif);
                mutations.push_back(mut);

            }

            free(reg);
        }

        LOGGER.info << "Total number of sequences: " << count << std::endl;

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
