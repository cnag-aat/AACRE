GENCCODE_ASSEMBLY_TYPE=-a gcc
SO=so
SOBJ=so
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.68.1
COMPILE=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/_build_env/bin/x86_64-conda-linux-gnu-cc -DNDEBUG -D_FORTIFY_SOURCE=2 -O2 -isystem /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include -D_REENTRANT  -DU_HAVE_ELF_H=1 -DU_HAVE_STRTOD_L=1 -DU_HAVE_XLOCALE_H=1 -DU_HAVE_STRING_VIEW=1  -DU_ATTRIBUTE_DEPRECATED= -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include -fdebug-prefix-map=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/work=/usr/local/src/conda/icu-68.1 -fdebug-prefix-map=/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -c
LIBFLAGS=-I/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include -DPIC -fPIC
GENLIB=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/_build_env/bin/x86_64-conda-linux-gnu-cc -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include -fdebug-prefix-map=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/work=/usr/local/src/conda/icu-68.1 -fdebug-prefix-map=/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,-rpath,/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib -Wl,-rpath-link,/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib -L/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib  -shared -Wl,-Bsymbolic
LDICUDTFLAGS=-nodefaultlibs -nostdlib
LD_SONAME=-Wl,-soname -Wl,
RPATH_FLAGS=
BIR_LDFLAGS=-Wl,-Bsymbolic
AR=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/_build_env/bin/x86_64-conda-linux-gnu-ar
ARFLAGS=r
RANLIB=/home/conda/feedstock_root/build_artifacts/icu_1605893451203/_build_env/bin/x86_64-conda-linux-gnu-ranlib
INSTALL_CMD=/usr/bin/install -c
