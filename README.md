samtools
========

This is *not* the official development repository for samtools.

This is our version of samtools, based on the original version, on which we have added an interactive BAM viewer tool.

The official repository is here:
- [samtools](https://github.com/samtools/samtools): mpileup and other tools for handling SAM, BAM, CRAM

To build this tool, you will need to checkout the htslib project as ../htslib

The htslib can be found here:
- [htslib](https://github.com/samtools/htslib): C-library for handling high-throughput sequencing data

See also http://github.com/samtools/

### Building Samtools

See [INSTALL](INSTALL) for complete details.
[Release tarballs][download] contain generated files that have not been
committed to this repository, so building the code from a Git repository
requires extra steps:

```sh
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
make install
```

By default, this will build against an HTSlib source tree in `../htslib`.
You can alter this to a source tree elsewhere or to a previously-installed
HTSlib by configuring with `--with-htslib=DIR`.

[download]: http://www.htslib.org/download/
