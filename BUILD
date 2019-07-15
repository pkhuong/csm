cc_library(
    name = "csm",
    srcs = ["csm.c"],
    hdrs = ["csm.h"],
    visibility = ["//visibility:public"],
    deps = [],
)

cc_test(
    name = "csm_test",
    srcs = ["csm.c", "csm.h"],
    copts = ["-DTEST_CSM"],
    deps = [],
)

cc_test(
    name = "csm-stat_test",
    srcs = ["csm-stat_test.cc"],
    size = "enormous",  # we need a ton of data points
    deps = [
        ":csm",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "csm-quantile-stat_test",
    srcs = ["csm-quantile-stat_test.cc"],
    size = "enormous",
    shard_count = 5,
    deps = [
        ":csm",
        "@com_google_googletest//:gtest_main",
    ],
)
