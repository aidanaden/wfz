const std = @import("std");

pub fn build(b: *std.Build) void {
    // Allow the user to enable or disable Tracy support with a build flag
    const tracy_enabled = b.option(
        bool,
        "tracy",
        "Build with Tracy support.",
    ) orelse true;

    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const zigimg_dep = b.dependency("zigimg", .{
        .target = target,
        .optimize = optimize,
    });

    // Get the Tracy dependency
    const tracy_dep = b.dependency("tracy", .{
        .target = target,
        .optimize = optimize,
        // ...
    });

    const exe_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    const exe = b.addExecutable(.{
        .name = "wfz",
        .root_module = exe_mod,
    });
    exe.root_module.addImport("zigimg", zigimg_dep.module("zigimg"));

    exe.root_module.addImport("tracy", tracy_dep.module("tracy"));
    // The user asked to enable Tracy, use the real implementation otherwise use dummy implementation
    exe.root_module.addImport("tracy_impl", tracy_dep.module(if (tracy_enabled) "tracy_impl_enabled" else "tracy_impl_disabled"));

    // Install
    b.installArtifact(exe);

    // Run
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // Tests
    const exe_unit_tests = b.addTest(.{
        .root_module = exe_mod,
    });
    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_exe_unit_tests.step);
}
