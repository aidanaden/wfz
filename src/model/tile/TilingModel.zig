const std = @import("std");
const Allocator = std.mem.Allocator;

const BaseModel = @import("../Base.zig");
const ImageProcessor = @import("../../ImageProcessor.zig");

const TileSet = @import("TileSet.zig");

const TilingModel = @This();

base: BaseModel,

tileset: TileSet,

/// All possible neighbor patterns for each cell,
/// values are initialised ONCE at the start then
/// remain untouched
propagator: [][]std.DynamicBitSet,

/// Processes and stores tile image pixel data
img_processor: ImageProcessor,

/// All possible tile pixel idxs
///
/// Note that this is DIFFERENT from `img_processor.img_pixel_idxs` since
/// that stores the original source tile image pixel idxs, but `tile_pixel_idxs`
/// contains pixel idxs of all possible permutations of source tile images
tile_pixel_idxs: [][]u32,

weight_log_weights: []f32,

tile_w: usize,
tile_h: usize,

allocator: Allocator,

pub fn init(
    params: *const BaseModel.InitParams,
    allocator: Allocator,
) !TilingModel {
    const img_dir_path = try std.fmt.allocPrintZ(allocator, "./src/assets/tilesets/{s}", .{params.name});
    std.debug.print("opening directory: {s}\n", .{img_dir_path});
    var img_dir = try std.fs.cwd().openDirZ(img_dir_path, .{
        .iterate = true,
    });
    defer img_dir.close();

    var img_paths = std.ArrayList([]const u8).init(allocator);
    defer img_paths.deinit();

    // Link image name to idx
    var img_name_idx = std.StringArrayHashMap(u32).init(allocator);
    defer img_name_idx.deinit();

    var iter = img_dir.iterate();
    while (try iter.next()) |e| {
        try img_name_idx.put(e.name, @as(u32, @intCast(img_paths.items.len)));

        const img_path = try std.fmt.allocPrintZ(allocator, "./src/assets/tilesets/{s}/{s}", .{ params.name, e.name });
        try img_paths.append(img_path);
    }

    var img_processor = try ImageProcessor.init(img_paths.items, allocator);
    defer img_processor.deinit();

    const tile_h = img_processor.images[0].height;
    const tile_w = img_processor.images[0].width;

    for (img_processor.images, 0..) |image, i| {
        const num_pixels = image.pixels.len();
        const img_path = img_paths.items[i];
        std.debug.print("image path: {s}, num pixels: {}\n", .{ img_path, num_pixels });
    }

    const tileset_path = try std.fmt.allocPrintZ(allocator, "./src/assets/tilesets/{s}.json", .{params.name});
    const tileset = try TileSet.fromJsonFile(tileset_path, allocator);

    try tileset.print_tile_neighbors();

    var weight_log_weights = try allocator.alloc(f32, tileset.source_tiles.len);

    var weight_total: f32 = 0;
    var weight_log_weight_total: f32 = 0.0;

    for (tileset.source_tiles, 0..) |tile, i| {
        weight_total += tile.weight;
        const weight_log_weight: f32 = @as(f32, @log(tile.weight));
        weight_log_weight_total += weight_log_weight;
        weight_log_weights[i] = weight_log_weight;
    }

    const base = try BaseModel.init(
        params,
        params.output_w / tile_w,
        params.output_h / tile_h,
        tileset.source_tiles.len,
        weight_total,
        weight_log_weight_total,
        allocator,
    );

    var propagator = try allocator.alloc([]std.DynamicBitSet, tileset.source_tiles.len);
    for (tileset.source_tiles, 0..) |tile, tile_i| {
        var neighbors = try allocator.alloc(std.DynamicBitSet, BaseModel.Direction.len());
        for (0..BaseModel.Direction.len()) |dir_i| {
            neighbors[dir_i] = tile.neighbors[dir_i];
        }
        propagator[tile_i] = neighbors;
    }

    return TilingModel{
        .tileset = tileset,
        .base = base,
        .weight_log_weights = weight_log_weights,
        .propagator = propagator,
        .img_processor = img_processor,
        .tile_h = tile_h,
        .tile_w = tile_w,
        // For now we ONLY use the original source tiles
        .tile_pixel_idxs = img_processor.img_pixel_idxs,
        .allocator = allocator,
    };
}

pub fn deinit(self: *TilingModel) void {
    // Propagator
    for (self.propagator) |neighbors| {
        // for (neighbors) |*n| {
        //     n.deinit();
        // }
        self.allocator.free(neighbors);
    }
    self.allocator.free(self.propagator);

    self.tileset.deinit();
    self.base.deinit();
    self.img_processor.deinit();

    self.allocator.free(self.weight_log_weights);
}

const tracy = @import("tracy");
pub fn propagate(self: *TilingModel) !void {
    std.debug.print("propagate\n", .{});
    const zone = tracy.Zone.begin(.{
        .name = "propagate",
        .src = @src(),
        .color = .orange,
        .callstack_depth = 10,
    });
    defer zone.end();

    while (self.base.propagation_stack.pop()) |idx| {
        self.base.on_stack.setValue(idx, false);

        const x: usize = @mod(idx, self.base.output_cells_w);
        const y: usize = @divFloor(idx, self.base.output_cells_w);

        // Possible patterns of cell at `idx`
        const source_patterns = self.base.possibilities[idx];

        for (std.meta.tags(BaseModel.Direction), 0..) |d, i| {
            // Calculate opp direction
            const opp_d = (i + 2) % 4;
            const dx, const dy = d.point();

            const neig_x: ?u32 = blk: {
                var neig_x: i32 = @as(i32, @intCast(x)) + dx;
                if (neig_x < 0) {
                    if (!self.base.periodic) {
                        break :blk null;
                    }
                    neig_x += @as(i32, @intCast(self.base.output_cells_w));
                }
                if (neig_x >= self.base.output_cells_w) {
                    if (!self.base.periodic) {
                        break :blk null;
                    }
                    neig_x -= @as(i32, @intCast(self.base.output_cells_w));
                }
                break :blk @as(u32, @intCast(neig_x));
            };

            const neig_y: ?u32 = blk: {
                var neig_y: i32 = @as(i32, @intCast(y)) + dy;
                if (neig_y < 0) {
                    if (!self.base.periodic) {
                        break :blk null;
                    }
                    neig_y += @as(i32, @intCast(self.base.output_cells_h));
                }
                if (neig_y >= self.base.output_cells_h) {
                    if (!self.base.periodic) {
                        break :blk null;
                    }
                    neig_y -= @as(i32, @intCast(self.base.output_cells_h));
                }
                break :blk @as(u32, @intCast(neig_y));
            };

            if (neig_x == null or neig_y == null) {
                continue;
            }

            const neig_idx = neig_x.? + neig_y.? * self.base.output_cells_w;
            const neig_patterns = self.base.possibilities[neig_idx];

            // For every possible pattern at neighbor idx
            var remove_idx = try std.ArrayList(usize).initCapacity(self.base.allocator, source_patterns.capacity());
            defer remove_idx.deinit();

            var neig_pattern_iter = neig_patterns.iterator(.{});
            while (neig_pattern_iter.next()) |neig_pattern_idx| {
                // For every possible neighbor pattern in current direction,
                // check if neighbor pattern is possible by checking if
                // neighbor pattern idx is present as a possible pattern
                // for cell at `idx`
                //
                // If not, remove pattern from possible patterns of cell at `idx`
                const required_neig_pattern_neigbors = self.propagator[neig_pattern_idx][opp_d];
                const is_supported = blk: {
                    var neigbors_iter = required_neig_pattern_neigbors.iterator(.{});
                    while (neigbors_iter.next()) |neigbor_idx| {
                        if (source_patterns.isSet(neigbor_idx)) {
                            break :blk true;
                        }
                    }
                    break :blk false;
                };

                // Remove neighbor pattern if current idx not supported
                if (!is_supported) {
                    remove_idx.appendAssumeCapacity(neig_pattern_idx);
                }
            }

            for (remove_idx.items) |remove_tile_idx| {
                try self.remove_tile(neig_idx, remove_tile_idx);
            }
        }
    }
}

/// Single iteration of collapsing probabilities
///
/// Returns completion status
pub fn collapse(
    self: *TilingModel,
    // Reusable allocated buffer to skip re-allocating each iteration
    random: std.Random,
) !BaseModel.Collapsed {
    var possible_tile_weights = try self.allocator.alloc(f32, self.tileset.source_tiles.len);
    defer self.allocator.free(possible_tile_weights);

    std.debug.print("collapse\n", .{});
    const zone = tracy.Zone.begin(.{
        .name = "collapse",
        .src = @src(),
        .color = .blue,
        .callstack_depth = 10,
    });
    defer zone.end();

    var min: f32 = 1000.0;
    var argminx: i32 = -1;
    var argminy: i32 = -1;

    // Find the point with minimum entropy (adding a little noise for randomness)
    for (0..self.base.output_cells_h) |y| {
        for (0..self.base.output_cells_w) |x| {
            const on_boundary = !self.base.periodic and (x > self.base.output_cells_w or y > self.base.output_cells_h);
            if (on_boundary) {
                continue;
            }

            const idx = x + y * self.base.output_cells_w;
            const remaining = self.base.possibilities_remaining[idx];

            // Contradiction: cell has 0 possible patterns, break
            if (remaining == 0) {
                std.debug.print("collapse: finished unsuccessfuly, contradiction found.\n", .{});
                return .contradiction;
            }

            if (remaining > 1) {
                const weights_sum = self.base.weight_sums[idx];
                const weight_log_weight_sum = self.base.weight_log_weight_sums[idx];
                const entropy = @as(f32, @log(weights_sum)) - (weight_log_weight_sum / weights_sum);

                // Add noise and check if this is the new minimum.
                // NOTE: lowest entropy = highest probability of a possible pattern
                const noise: f32 = 0.000001 * random.float(f32);
                if (entropy > 0 and (entropy + noise) < min) {
                    min = entropy + noise;
                    argminx = @as(i32, @intCast(x));
                    argminy = @as(i32, @intCast(y));
                }
            }

            if (remaining < 1) {
                std.debug.print("collapse: remaining not > 1 or not 0: {}\n", .{remaining});
            }
        }
    }

    // Check if the generation is complete.
    if (argminx == -1 and argminy == -1) {
        std.debug.print("collapse: finished successfully\n", .{});
        return .success; // Finished successfully
    }

    const argmin_idx = @as(usize, @intCast(argminx + argminy * @as(i32, @intCast(self.base.output_cells_w))));

    // Assign probability of all patterns to the probability
    // of the pattern w the lowest entropy
    var argmin_possible_iter = self.base.possibilities[argmin_idx].iterator(.{});
    while (argmin_possible_iter.next()) |tile_idx| {
        const tile = self.tileset.source_tiles[tile_idx];
        possible_tile_weights[tile_idx] = tile.weight;
    }

    const random_idx = try random_tile_idx(possible_tile_weights, random.float(f32));

    var remove = try std.ArrayList(usize).initCapacity(self.allocator, possible_tile_weights.len);
    defer remove.deinit();
    var iter = self.base.possibilities[argmin_idx].iterator(.{});
    while (iter.next()) |pattern_idx| {
        if (pattern_idx != random_idx) {
            remove.appendAssumeCapacity(pattern_idx);
        }
    }

    std.debug.print("removing {} patterns from argmin_idx {}\n", .{ remove.items.len, argmin_idx });
    for (remove.items) |idx| {
        try self.remove_tile(argmin_idx, idx);
    }
    return .failed;
}

/// Return the first index `x` where `sum of values [0 .. x - 1] > r`.
///
/// This is a form of doing weighted RNG, where the elements in `[0 .. x - 1]`
/// with higher value have a higher chance of being chosen, since greater
/// values have higher chance of being greater than `r`
fn random_tile_idx(tile_weights: []f32, r: f32) !usize {
    const sum = blk: {
        var s: f32 = 0.0;
        for (tile_weights) |weight| {
            s += weight;
        }
        if (s == 0) {
            @memset(tile_weights, 1.0);
            s = @as(f32, @floatFromInt(tile_weights.len));
        }
        break :blk s;
    };

    var x: f32 = 0.0;
    var last_key: usize = undefined; // Keep track of the last valid key
    for (tile_weights, 0..) |*weight, i| {
        last_key = i;
        const normalised: f32 = weight.* / sum;
        weight.* = normalised;
        x += normalised;
        if (r <= x) {
            return i;
        }
    }
    return last_key;
}

pub fn remove_tile(self: *TilingModel, idx: usize, tile_idx: usize) !void {
    // Already removed, break
    if (!self.base.possibilities[idx].isSet(tile_idx)) {
        return;
    }

    self.base.possibilities[idx].unset(tile_idx);
    self.base.possibilities_remaining[idx] -= 1;

    const tile = self.tileset.source_tiles[tile_idx];
    self.base.weight_sums[idx] -= tile.weight;

    const weight_log_weight = self.weight_log_weights[tile_idx];
    self.base.weight_log_weight_sums[idx] -= weight_log_weight;

    const added_to_stack = self.base.on_stack.isSet(idx);
    if (!added_to_stack) {
        try self.base.propagation_stack.append(idx);
        self.base.on_stack.set(idx);
    }
}

pub fn clear(self: *TilingModel) void {
    var weight_total: f32 = 0.0;
    var weight_log_weight_total: f32 = 0.0;
    for (self.tileset.source_tiles, 0..) |tile, tile_idx| {
        weight_total += tile.weight;
        // Note: self.weight_log_weights is static and already calculated.
        // We can just sum its values.
        weight_log_weight_total += self.weight_log_weights[tile_idx];
    }
    const num_patterns = self.tileset.source_tiles.len;
    self.base.clear(
        num_patterns,
        weight_total,
        weight_log_weight_total,
    );
}

// pub fn collapse_ground(self: *TilingModel) !void {
// if (self.ground_patterns.keys().len == 0) {
//     return;
// }
//
// const ground_y: usize = self.base.output_cells_h - 1;
// var remove = try std.ArrayList(usize).initCapacity(
//     self.allocator,
//     self.pattern_weights.keys().len,
// );
// defer remove.deinit();
//
// for (0..self.base.output_cells_w) |x| {
//     const ground_idx: usize = x + ground_y * self.base.output_cells_w;
//     const ground_possibilities = self.base.possibilities[ground_idx];
//
//     var ground_possiblities_iter = ground_possibilities.iterator(.{});
//     while (ground_possiblities_iter.next()) |pattern_idx| {
//         const pattern_id = self.idx_to_pattern[pattern_idx];
//         if (self.ground_patterns.get(pattern_id) == null) {
//             remove.appendAssumeCapacity(pattern_idx);
//         }
//     }
//
//     for (remove.items) |remove_idx| {
//         try self.remove_pattern(ground_idx, remove_idx);
//     }
//     remove.clearRetainingCapacity();
// }
// }

// Given `idx`, generates output pixels for the next `cell_h` x `cell_w` grid of pixels
// pub fn render_cell(self: *TilingModel, cell_idx: usize, output_pixel_ids: []usize) !void {
//     const cell_patterns = self.base.possibilities[cell_idx];
//     const pattern_idx = cell_patterns.findFirstSet() orelse {
//         return error.NoPossiblePatterns;
//     };
// }
//

/// Generates output pixels for the cell at `cell_idx`
pub fn render_cell(self: *TilingModel, cell_idx: usize, output_pixel_idxs: []u32) !void {
    const cell_tiles = self.base.possibilities[cell_idx];
    const tile_idx = cell_tiles.findFirstSet() orelse {
        return error.NoPossiblePatterns;
    };

    const pixel_idxs = self.img_processor.img_pixel_idxs[tile_idx];

    std.debug.assert(pixel_idxs.len == self.tile_h * self.tile_w);
    std.debug.assert(output_pixel_idxs.len > pixel_idxs.len);

    // NOTE: output pixel is the FIRST (top left) pixel of the pattern
    // since for overlap model cell_h and cell_w is 1x1
    for (0..self.tile_h) |dy| {
        for (0..self.tile_w) |dx| {
            const idx = dx + dy * self.tile_w;
            const pixel_idx: u32 = pixel_idxs[idx];
            const output_idx = cell_idx + idx;
            output_pixel_idxs[output_idx] = pixel_idx;
        }
    }
}

/// Output slice MUST be deinitialised
pub fn render_output(self: *TilingModel) ![]u32 {
    const output_pixel_idxs: []u32 = try self.allocator.alloc(u32, self.base.output_cells_w * self.base.output_cells_h);
    for (0..self.base.output_cells_h) |y| {
        for (0..self.base.output_cells_w) |x| {
            const cell_idx = x + self.base.output_cells_w * y;
            try self.render_cell(cell_idx, output_pixel_idxs);
        }
    }
    return output_pixel_idxs;
}
