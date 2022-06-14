./misc/condorJobMaker.py \
        cmsConfigs/data2018Parking_v1.py \
        fileList/Run2018D_ParkingBPH3_AOD_20Jun2021_UL2018-v1_270001.fls \
        results/data/Run2018D_ParkingBPH3/part2p3 \
        10000 \
        1 \
        -1 \
        dataBPHParking3D \
        100

./misc/condorJobMaker.py \
        cmsConfigs/mcUL2018_v1.py \
        fileList/BuToJpsiK_aodsim.fls \
        results/MC/BuToJpsiK/part2p3 \
        200 \
        1 \
        -1 \
        buToJpsiK \
        50


./misc/condorJobMaker.py \
        cmsConfigs/mcUL2018_v1.py \
        fileList/bsTommg_UL2018_RECO.fls \
        results/MC/BsToMuMuGamma/part2p3 \
        200 \
        1 \
        -1 \
        bsToMuMuGamma \
        75

./misc/condorJobMaker.py \
        cmsConfigs/mcUL2018_v1.py \
        fileList/BsToJPsiGamma.fls \
        results/MC/BsToJPsiGamma/part2p3 \
        200 \
        1 \
        -1 \
        bsToJPsiGamma \
        75
