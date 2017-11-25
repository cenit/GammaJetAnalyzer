
PKGTOOLS_TAG=V00-30-XX;

PR_TESTS=1;
PROD_ARCH=1;
ENABLE_DEBUG=1;
NO_IB=1;

#RELEASE_QUEUE=CMSSW_8_4_X;
#CMSDIST_TAG=IB/CMSSW_8_4_X/gcc530;
#CMSSW_X_Y_Z=CMSSW_8_4_X
RELEASE_QUEUE=CMSSW_9_0_X;
CMSDIST_TAG=IB/CMSSW_9_0_X/gcc530;
CMSSW_X_Y_Z=CMSSW_9_0_X
SCRAM_ARCH=slc6_amd64_gcc530;
ARCH=slc6_amd64_gcc530

mkdir -p $CMSSW_X_Y_Z-build
cd $CMSSW_X_Y_Z-build
QUEUE=$(echo $CMSSW_X_Y_Z | sed -e 's/\(CMSSW_[0-9][0-9]*_[0-9][0-9]*\).*/\1_X/')
eval $(curl -k -s https://raw.githubusercontent.com/cms-sw/cms-bot/master/config.map | grep "SCRAM_ARCH=$ARCH;" | grep "RELEASE_QUEUE=$QUEUE;")
git clone -b $CMSDIST_TAG git@github.com:cms-sw/cmsdist.git CMSDIST
git clone -b $PKGTOOLS_TAG git@github.com:cms-sw/pkgtools.git PKGTOOLS

sh -e PKGTOOLS/scripts/prepare-cmsdist $CMSSW_X_Y_Z $ARCH

time PKGTOOLS/cmsBuild --architecture=$ARCH --builders 4 -j 16 build cms-git-tools
time PKGTOOLS/cmsBuild --architecture=$ARCH --builders 4 -j 16 build cmssw

