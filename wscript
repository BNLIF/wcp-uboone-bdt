#!/usr/bin/env python


TOP = '.'
APPNAME = 'Junk'

def options(opt):
    opt.load("wcb")

def configure(cfg):
    cfg.load("wcb")

    if cfg.env["HAVE_VLNEVAL"]:
        cfg.check_boost(lib = 'program_options')

    # boost 1.59 uses auto_ptr and GCC 5 deprecates it vociferously.
    cfg.env.CXXFLAGS += ['-Wno-deprecated-declarations']
    cfg.env.CXXFLAGS += ['-Wall', '-Wno-unused-local-typedefs', '-Wno-unused-function']
    # cfg.env.CXXFLAGS += ['-Wpedantic', '-Werror']


def build(bld):
    import os
    bld.load('wcb')

    use = [ 'ROOTSYS' ]

    if bld.env["HAVE_VLNEVAL"]:
        use += [ 'VLNEVAL', 'BOOST' ]

    bld.smplpkg('WCPuBooNE_BDT_APP', use = use)

    # Set test-runner LD_LIBRARY_PATH so tests use the freshly-built library
    # (RUNPATH < LD_LIBRARY_PATH, so build/ must come first to beat bdt_install).
    build_dir = bld.path.get_bld().abspath()
    root_libs = ':'.join(bld.env.LIBPATH_ROOTSYS or [])
    test_env = dict(os.environ)
    test_env['LD_LIBRARY_PATH'] = (
        build_dir + ':' + root_libs + ':' + test_env.get('LD_LIBRARY_PATH', '')
    )
    for g in bld.groups:
        for tg in g:
            if getattr(tg, 'ut_cwd', None) is not None:
                tg.ut_env = test_env

