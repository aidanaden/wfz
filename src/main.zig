const std = @import("std");
const Allocator = std.mem.Allocator;
const zigimg = @import("zigimg");
const Pattern = @import("pattern.zig").Pattern;
const SourcePattern = Pattern(u32);

const OverlapModel = struct {
    const Options = struct {
        N: u32,
        periodic_input: bool,
        periodic_output: bool,
        symmetry: u32,
        ground: u32,
        const DEFAULTS: Options = .{
            .N = 2,
            .periodic_input = true,
            .periodic_output = false,
            .symmetry = 8,
            .ground = 0,
        };
    };

    /// Output image width
    width: u32,
    /// Output image height
    height: u32,
    /// Pattern size (NxN)
    N: u32,
    periodic: bool,
    symmetry: u32,
    patterns: []SourcePattern,

    /// Id of the specific pattern to use as the bottom of the
    /// generation, a value of 0 means that this is unset
    ground: u32,

    const Self = @This();
    fn init(
        width: u32,
        height: u32,
        name: []const u8,
        options: ?Options,
    ) Self {
        const _options = options orelse Options.DEFAULTS;
        return Self{
            .width = width,
            .height = height,
            .name = name,
            .N = _options.N,
            .periodic = _options.periodic_output,
            .symmetry = _options.symmetry,
            .ground = _options.ground,
        };
    }
};

pub fn main() !void {
    // Init allocator
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    const allocator = arena.allocator();

    const img_path = "./src/assets/City.png";
    var image = try zigimg.Image.fromFilePath(allocator, img_path);
    defer image.deinit();

    const pixels = image.pixels.rgba32;
    // TODO: dynamically set this based on return type of `image.pixels`
    const RgbType = zigimg.color.Rgba32;

    // Mapping of color_id -> rgb value
    var id_to_pixel = std.AutoArrayHashMap(u32, RgbType).init(allocator);
    defer id_to_pixel.deinit();

    // Mapping of rgb value -> color_id
    var pixel_to_id = std.AutoArrayHashMap(RgbType, u32).init(allocator);
    defer pixel_to_id.deinit();

    // Initial starting color id value
    var color_id: u32 = 1;

    // All pixels of the source image encoded in color_id
    var pixel_ids = try std.ArrayList(u32).initCapacity(allocator, pixels.len);
    defer pixel_ids.deinit();

    // Encode every pixel into a color id
    for (pixels) |pixel| {
        const id = pixel_to_id.get(pixel);
        if (id == null) {
            try pixel_to_id.put(pixel, color_id);
            try id_to_pixel.put(color_id, pixel);
            try pixel_ids.append(color_id);
            color_id += 1;
        } else {
            try pixel_ids.append(id.?);
        }
    }

    const num_colors: u32 = @intCast(pixel_to_id.keys().len);
    const pattern_encoder = SourcePattern.BaseEncoder.init(3, num_colors);
    var patterns = std.ArrayList(SourcePattern).init(allocator);
    defer {
        for (patterns.items) |pattern| {
            pattern.deinit();
        }
        patterns.deinit();
    }

    // Generate patterns based on NxN tiling window
    for (0..image.height) |y| {
        for (0..image.width) |x| {
            const pattern = try SourcePattern.init(
                x,
                y,
                3,
                pixel_ids.items,
                image.height,
                image.width,
                allocator,
            );
            const patternId = pattern_encoder.encode(pattern.values.items);
            std.debug.print("pattern ({}, {}): id: {} \n", .{ x, y, patternId });
            for (pattern.values.items, 0..) |v, i| {
                const pixel = id_to_pixel.get(v);
                if (pixel == null) {
                    std.debug.print("i: {}, pixel_id: {}: pixel_value: UNKNOWN\n", .{ i, v });
                } else {
                    std.debug.print("i: {}, pixel_id: {}: pixel_value: {}\n", .{ i, v, pixel.? });
                }
            }
            try patterns.append(pattern);
        }
    }

    std.debug.print("num_colors: {}, width: {}, height: {}, pattern count: {}\n", .{ num_colors, image.width, image.height, patterns.items.len });
}
