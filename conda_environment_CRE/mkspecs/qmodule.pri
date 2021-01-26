QMAKE_CFLAGS_WARN_ON += -Wno-expansion-to-defined
QMAKE_CXXFLAGS_WARN_ON += -Wno-expansion-to-defined
EXTRA_DEFINES += _X_INLINE=inline XK_dead_currency=0xfe6f _FORTIFY_SOURCE=2 XK_ISO_Level5_Lock=0xfe13 FC_WEIGHT_EXTRABLACK=215 FC_WEIGHT_ULTRABLACK=FC_WEIGHT_EXTRABLACK GLX_GLXEXT_PROTOTYPES
EXTRA_INCLUDEPATH += /projects/66d93023-00f0-4c12-8a25-5d6d4e486740/miniforge3/conda-bld/qt_1608618950355/work/openssl_hack/include /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include
EXTRA_LIBDIR += /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/x86_64-conda-linux-gnu/sysroot/usr/lib64
!host_build|!cross_compile {
    QMAKE_LFLAGS+=-Wl,-rpath,/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib -Wl,-rpath-link,/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib -L/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib
}
QT_CPU_FEATURES.x86_64 = mmx sse sse2
QT.global_private.enabled_features = sse2 alloca_h alloca dbus dbus-linked gui network posix_fallocate reduce_exports reduce_relocations sql system-zlib testlib widgets xml
QT.global_private.disabled_features = alloca_malloc_h android-style-assets avx2 private_tests gc_binaries libudev release_tools stack-protector-strong
PKG_CONFIG_EXECUTABLE = /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/bin/pkg-config
QMAKE_LIBS_DBUS = /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib/libdbus-1.so
QMAKE_INCDIR_DBUS = /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include/dbus-1.0 /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib/dbus-1.0/include
QT_COORD_TYPE = double
QMAKE_LIBS_ZLIB = /home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib/libz.so
CONFIG += use_gold_linker sse2 aesni compile_examples enable_new_dtags largefile optimize_size precompile_header rdrnd shani sse3 ssse3 sse4_1 sse4_2 x86SimdAlways
QT_BUILD_PARTS += libs tools
QT_HOST_CFLAGS_DBUS += -I/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/include/dbus-1.0 -I/home/devel/jgomez/AACRE_pipeline.V01/conda_environment_CRE/lib/dbus-1.0/include
