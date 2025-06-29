/// tracy imports
pub const tracy = @import("tracy");
pub const tracy_impl = @import("tracy_impl");
pub const tracy_options: tracy.Options = .{
    .default_callstack_depth = 8,
    .verbose = true,
};

const std = @import("std");
const Allocator = std.mem.Allocator;
const BaseModel = @import("./model/Base.zig");
const OverlapModel = @import("./model/overlap/OverlapModel.zig");
const ImageProcessor = @import("ImageProcessor.zig");

const SourcePattern = @import("./model/overlap/Pattern.zig").Pattern(u32);
const PatternId = SourcePattern.BaseEncoder.EncodedId;

const TileSet = @import("./model/tile/TileSet.zig");
const TilingModel = @import("./model/tile/TilingModel.zig");

pub fn main() !void {
    const zone = tracy.Zone.begin(.{
        .name = "Main zone",
        .src = @src(),
        .color = .tomato,
    });
    defer zone.end();

    // var tracy_allocator: tracy.Allocator = .{ .parent = std.heap.page_allocator };

    // Init allocator
    // var debug = std.heap.DebugAllocator(.{ .verbose_log = true }).init;
    // defer {
    //     const leaks = debug.deinit();
    //     std.log.debug("leaks: {any}\n", .{leaks});
    // }
    // const debug_allocator = debug.allocator();
    // var arena = std.heap.ArenaAllocator.init(debug_allocator);

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const name = "castle";
    // const name = "bomberman-stage-3-08";
    // const name = "3Bricks";
    // const name = "city";
    const img_path = try std.fmt.allocPrintZ(allocator, "./src/assets/samples/{s}.png", .{name});
    const output_fname = try std.fmt.allocPrint(allocator, "{s}_output.png", .{name});

    // Set test parameters
    const params = BaseModel.InitParams{
        .N = 3,
        .periodic = true,
        .output_w = 256, // 340,
        .output_h = 256, // 512,
        .symmetry = 8,
        .ground = false,
        .img_path = img_path,
        .output_fname = output_fname,
        .name = name,
    };
    // load asset(s)
    //
    // build relation probabilities
    //
    // generate

    // tracy.message(.{
    //     .text = "overlapInit",
    //     .color = .light_green,
    // });

    var tiling = try TilingModel.init(&params, allocator);
    defer tiling.deinit();

    var prng = std.Random.DefaultPrng.init(@as(u64, @intCast(std.time.nanoTimestamp())));
    var random = prng.random();

    var collapsed: BaseModel.Collapsed = .failed;
    while (collapsed != .success) {
        tiling.clear();
        random = prng.random();
        collapsed = try generate(&tiling, random);
    }

    const output_pixel_idxs = try tiling.render_output();
    defer allocator.free(output_pixel_idxs);

    var output_img = try tiling.img_processor.generate(
        params.output_w,
        params.output_h,
        output_pixel_idxs,
    );
    defer output_img.deinit();

    std.debug.print("image pixels: {any}\n", .{output_img.imageByteSize()});
    const output_fpath = try std.fmt.allocPrintZ(allocator, "./src/assets/tiling/{s}", .{params.output_fname});
    try output_img.writeToFilePath(output_fpath, .{ .png = .{} });

    // var overlap = try OverlapModel.init(
    //     &params,
    //     allocator,
    // );
    // var collapsed: OverlapModel.Collapsed = .failed;
    // while (collapsed != .success) {
    //     overlap.clear();
    //     std.debug.print("Unique patterns found: {}\n", .{overlap.pattern_weights.keys().len});
    //     random = prng.random();
    //     collapsed = try generate(&overlap, random);
    // }
    //
    // const output_pixel_idxs = try overlap.render_output();
    // defer allocator.free(output_pixel_idxs);
    //
    // var output_img = try overlap.img_processor.generate(
    //     params.output_w,
    //     params.output_h,
    //     output_pixel_idxs,
    // );
    // defer output_img.deinit();
    //
    // std.debug.print("image pixels: {any}\n", .{output_img.imageByteSize()});
    // const output_fpath = try std.fmt.allocPrintZ(allocator, "./src/assets/samples/{s}", .{params.output_fname});
    // try output_img.writeToFilePath(output_fpath, .{ .png = .{} });
}

// This function uses compile-time duck typing. It will compile for any `model`
// that has a `base: BaseModel` field and a `propagate()` method.
pub fn generate(model: anytype, random: std.Random) !BaseModel.Collapsed {
    std.debug.print("generate\n", .{});

    // if (model.ground_patterns.keys().len > 0) {
    //     try model.collapse_ground();
    // }
    try model.propagate();

    var i: usize = 0;
    var collapsed: BaseModel.Collapsed = .failed;
    while (collapsed == .failed) {
        std.debug.print("iterate {}\n", .{i});

        collapsed = try model.collapse(random);
        try model.propagate();

        if (collapsed == .success) {
            return .success;
        }
        i += 1;
    }

    std.debug.print("collapsed: {}\n", .{collapsed});
    return collapsed;
}

test {
    _ = @import("./model/tile/TileSet.zig");
}
