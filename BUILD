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