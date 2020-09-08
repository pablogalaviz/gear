//
// Created by pablo on 6/8/20.
//

#ifndef GEAR_TELOMEREREGION_H
#define GEAR_TELOMEREREGION_H

#include <genomeRegion.h>

namespace cmri {

    struct telomere_signature_t {
        unsigned int count;
        std::map<unsigned int,unsigned int> histogram;
        std::map<std::string, unsigned int> context;

        std::string serialize() const;
        void deserialize(const boost::property_tree::ptree &tree);

    };

    inline std::ostream& operator<<(std::ostream& result, const telomere_signature_t& rhs)
    {
        result << rhs.serialize();
        return result;
    }


    struct telomereRegion : public genomeRegion {

            unsigned int total_bases=0;

            std::map<std::string,telomere_signature_t> signature;

            std::string serialize() const override;
            void deserialize(const boost::property_tree::ptree &tree) override;


            virtual bool operator==(const telomereRegion &rhs) const {
                return genomeRegion::operator==(rhs) &&
                       total_bases == rhs.total_bases
                        ;
            }


        };

        inline std::ostream& operator<<(std::ostream& result, const telomereRegion& rhs)
        {
            result << rhs.serialize();
            return result;
        }


}


#endif //GEAR_TELOMEREREGION_H
