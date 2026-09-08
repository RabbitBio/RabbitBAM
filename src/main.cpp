#include <cstdio>
#include <cstring>

#include "rb_bam_cli.h"

static void usage(FILE *stream) {
    std::fprintf(stream,
        "RabbitBAM (sortedbam): BAM I/O modules from RabbitBin\n\n"
        "Usage: rabbitbam <command> [options]\n\n"
        "Commands:\n"
        "  sortbam  Coordinate-sort a BAM and optionally write its BAI index\n"
        "  bai      Build a BAI index for a coordinate-sorted BAM\n\n"
        "Run rabbitbam <command> --help for command options.\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage(stderr);
        return 1;
    }
    if (std::strcmp(argv[1], "--help") == 0 ||
        std::strcmp(argv[1], "-h") == 0) {
        usage(stdout);
        return 0;
    }
    if (std::strcmp(argv[1], "sortbam") == 0)
        return rb_cmd_sortbam(argc - 1, argv + 1);
    if (std::strcmp(argv[1], "bai") == 0)
        return rb_cmd_bai(argc - 1, argv + 1);

    std::fprintf(stderr, "Unknown command: %s\n\n", argv[1]);
    usage(stderr);
    return 1;
}
