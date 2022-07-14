#!/bin/sh

cases=('koss_highnox_bnzdiol_UNS-ALK' 'koss_highnox_bnzdiol_UNS-AR1' 'koss_highnox_bnzdiol_UNS-AR2' 'koss_highnox_meffral_UNS-ALK' 'koss_highnox_meffral_UNS-AR1' 'koss_highnox_meffral_UNS-AR2')

for cs in ${cases[@]}
do
    cp ./cases/saprc11b.rxn.${cs} saprc11b.rxn
    echo "saprc14_rev1" | ./prep.exe
    mv saprc14_rev1.ftn ../../../src/saprc_files/saprc14_rev1.f.${cs}
    mv saprc14_rev1.mod ../../../src/saprc_files/saprc14_rev1.mod.${cs}
#    mv saprc14_rev1.mod ./cases/saprc14_rev1.mod.${cs}
#    cp ./cases/saprc14_rev1.f.${cs} /data/n000/shantanu/src/airshed/ca24km/som_poa+ivoc/${cs}_24kmnest/saprc14_rev1.f
#    cp ./cases/saprc14_rev1.mod.${cs} /data/n000/shantanu/src/airshed/ca24km/som_poa+ivoc/${cs}_24kmnest/saprc14_rev1.mod
done
