const tracy = @import("tracy");

const std = @import("std");
const Allocator = std.mem.Allocator;
const Random = std.Random;

const SourcePattern = @import("./Pattern.zig").Pattern(u32);
const PatternId = SourcePattern.BaseEncoder.EncodedId;
const BaseModel = @import("../Base.zig");

const ImageProcessor = @import("../../ImageProcessor.zig");

const OverlapModel = @This();

const ImageInfo = struct {
    height_px: usize,
    width_px: usize,
    pixel_idxs: []u32,
};

base: BaseModel,

/// --- UNIQUE DATA ---
/// Data specific to the Overlapping model.
N: usize,

/// All possible neighbor patterns for each cell,
/// values are initialised ONCE at the start then
/// remain untouched
propagator: [][BaseModel.Direction.len()]std.DynamicBitSet,

/// Pattern related data
ground_patterns: std.AutoArrayHashMap(PatternId, void),
pattern_weights: std.AutoArrayHashMap(PatternId, u32),
weight_log_weights: std.AutoArrayHashMap(PatternId, f32),
pattern_to_idx: std.AutoArrayHashMap(PatternId, u32),
idx_to_pattern: []PatternId,

pattern_encoder: SourcePattern.BaseEncoder,

img_processor: ImageProcessor,

allocator: Allocator,

pub fn init(params: *const BaseModel.InitParams, allocator: Allocator) !OverlapModel {
    // Load source image
    var img_paths = [_][]const u8{params.img_path};
    const img_processor = try ImageProcessor.init(&img_paths, allocator);

    // Initialize pattern encoder (used to encode patterns from source image)
    const pattern_encoder = SourcePattern.BaseEncoder.init(
        params.N,
        img_processor.num_uniq_pixels,
    );

    // Construct patterns (encoded as pattern id) and weights based on input image
    const pattern_weights, const ground_patterns = try build_patterns(
        params,
        .{
            .height_px = img_processor.images[0].height,
            .width_px = img_processor.images[0].width,
            .pixel_idxs = img_processor.img_pixel_idxs[0],
        },
        &pattern_encoder,
        allocator,
    );

    var pattern_to_idx = std.AutoArrayHashMap(PatternId, u32).init(allocator);
    try pattern_to_idx.ensureTotalCapacity(pattern_weights.keys().len);

    var idx_to_pattern = try allocator.alloc(PatternId, pattern_weights.keys().len);
    var pattern_iter = pattern_weights.iterator();
    var dense_idx: u32 = 0;
    while (pattern_iter.next()) |entry| {
        const pattern_id = entry.key_ptr.*;
        // Map the sparse ID to the next available dense index
        pattern_to_idx.putAssumeCapacity(pattern_id, dense_idx);
        // Map the dense index back to the sparse ID
        idx_to_pattern[dense_idx] = pattern_id;
        dense_idx += 1;
    }

    var weight_log_weights = std.AutoArrayHashMap(PatternId, f32).init(allocator);
    try weight_log_weights.ensureTotalCapacity(pattern_weights.keys().len);

    var weight_total: u32 = 0;
    var weight_log_weight_total: f32 = 0.0;

    var pattern_weight_iter = pattern_weights.iterator();
    while (pattern_weight_iter.next()) |entry| {
        const count = entry.value_ptr.*;
        weight_total += count;
        const weight_log_weight: f32 = @as(f32, @floatFromInt(count)) * @as(f32, @log(@as(f32, @floatFromInt(count))));
        weight_log_weight_total += weight_log_weight;
        try weight_log_weights.put(entry.key_ptr.*, weight_log_weight);
    }

    const base = try BaseModel.init(
        params,
        params.output_w,
        params.output_h,
        pattern_weights.keys().len,
        weight_total,
        weight_log_weight_total,
        allocator,
    );

    const propagator = try init_propagator(
        pattern_weights,
        pattern_to_idx,
        &pattern_encoder,
        allocator,
    );

    return OverlapModel{
        .base = base,
        .N = params.N,
        .propagator = propagator,
        .allocator = allocator,
        .weight_log_weights = weight_log_weights,
        .pattern_weights = pattern_weights,
        .ground_patterns = ground_patterns,
        .idx_to_pattern = idx_to_pattern,
        .pattern_to_idx = pattern_to_idx,
        .pattern_encoder = pattern_encoder,
        .img_processor = img_processor,
    };
}

pub fn deinit(self: *OverlapModel) void {
    // Propagator
    for (self.propagator) |*neighbors| {
        for (neighbors) |*n| {
            n.deinit();
        }
    }
    self.allocator.free(self.propagator);

    // Pattern
    self.pattern_weights.deinit();
    self.allocator.free(self.idx_to_pattern);
    self.pattern_to_idx.deinit();
    self.ground_patterns.deinit();
    self.weight_log_weights.deinit();

    self.base.deinit();
    self.img_processor.deinit();
}

fn build_patterns(
    params: *const BaseModel.InitParams,
    image_info: ImageInfo,
    pattern_encoder: *const SourcePattern.BaseEncoder,
    allocator: Allocator,
) !struct {
    std.AutoArrayHashMap(PatternId, u32),
    std.AutoArrayHashMap(PatternId, void),
} {
    var pattern_weights = std.AutoArrayHashMap(PatternId, u32).init(allocator);
    var ground_patterns = std.AutoArrayHashMap(PatternId, void).init(allocator);

    const maxY, const maxW = if (params.periodic) .{
        image_info.height_px,
        image_info.width_px,
    } else .{
        image_info.height_px - params.N + 1,
        image_info.width_px - params.N + 1,
    };

    var pattern_ids = std.ArrayList(PatternId).init(allocator);
    defer pattern_ids.deinit();

    // Generate patterns based on NxN tiling window
    for (0..maxY) |y| {
        for (0..maxW) |x| {
            defer pattern_ids.clearRetainingCapacity();

            const pattern = try SourcePattern.init(
                x,
                y,
                params.N,
                image_info.pixel_idxs,
                image_info.height_px,
                image_info.width_px,
                allocator,
            );
            defer pattern.deinit();
            const original_id = pattern_encoder.encode(pattern.values.items);
            try pattern_ids.append(original_id);

            const reflect = try pattern.transform(.Reflect);
            defer reflect.deinit();
            const reflect_id = pattern_encoder.encode(reflect.values.items);
            try pattern_ids.append(reflect_id);

            const rotate = try pattern.transform(.Rotate);
            defer rotate.deinit();
            const rotate_id = pattern_encoder.encode(rotate.values.items);
            try pattern_ids.append(rotate_id);

            const reflect_rotate = try rotate.transform(.Reflect);
            defer reflect_rotate.deinit();
            const reflect_rotate_id = pattern_encoder.encode(reflect_rotate.values.items);
            try pattern_ids.append(reflect_rotate_id);

            const rotate_rotate = try rotate.transform(.Rotate);
            defer rotate_rotate.deinit();
            const rotate_rotate_id = pattern_encoder.encode(rotate_rotate.values.items);
            try pattern_ids.append(rotate_rotate_id);

            const reflect_rotate_rotate = try rotate_rotate.transform(.Reflect);
            defer reflect_rotate_rotate.deinit();
            const reflect_rotate_rotate_id = pattern_encoder.encode(reflect_rotate_rotate.values.items);
            try pattern_ids.append(reflect_rotate_rotate_id);

            const rotate_rotate_rotate = try rotate_rotate.transform(.Rotate);
            defer rotate_rotate_rotate.deinit();
            const rotate_rotate_rotate_id = pattern_encoder.encode(rotate_rotate_rotate.values.items);
            try pattern_ids.append(rotate_rotate_rotate_id);

            const reflect_rotate_rotate_rotate = try rotate_rotate_rotate.transform(.Reflect);
            defer reflect_rotate_rotate_rotate.deinit();
            const reflect_rotate_rotate_rotate_id = pattern_encoder.encode(reflect_rotate_rotate_rotate.values.items);
            try pattern_ids.append(reflect_rotate_rotate_rotate_id);

            for (0..params.symmetry) |i| {
                const pattern_id = pattern_ids.items[i];
                const decoded_pattern = try pattern_encoder.decode(pattern_id, allocator);
                defer decoded_pattern.deinit();

                // TESTING: display pattern rgb values
                // var cols: [9]RgbType = undefined; for (decoded_pattern.values.items, 0..) |v, i| {
                //     const col = id_to_pixel.get(v);
                //     if (col != null) {
                //         cols[i] = col.?;
                //     }
                // }

                // Increment pattern count
                const weight = pattern_weights.get(pattern_id) orelse 0;
                try pattern_weights.put(pattern_id, weight + 1);

                if (params.ground and (y == maxY - 1)) {
                    try ground_patterns.put(pattern_id, {});
                }
            }
        }
    }

    return .{ pattern_weights, ground_patterns };
}

const NumDirections = BaseModel.Direction.len();

fn init_propagator(
    pattern_weights: std.AutoArrayHashMap(PatternId, u32),
    pattern_to_idx: std.AutoArrayHashMap(PatternId, u32),
    pattern_encoder: *const SourcePattern.BaseEncoder,
    allocator: Allocator,
) ![][NumDirections]std.DynamicBitSet {
    var propagator = try allocator.alloc([NumDirections]std.DynamicBitSet, pattern_weights.keys().len);

    for (pattern_weights.keys()) |pattern_id1| {
        const pattern1 = try pattern_encoder.decode(pattern_id1, allocator);
        defer pattern1.deinit();
        var pattern_neigbors: [NumDirections]std.DynamicBitSet = undefined;

        inline for (std.meta.tags(BaseModel.Direction), 0..) |d, i| {
            const dx, const dy = d.point();
            var neighbors = try std.DynamicBitSet.initEmpty(allocator, pattern_weights.keys().len);
            for (pattern_weights.keys()) |pattern_id2| {
                var pattern2 = try pattern_encoder.decode(pattern_id2, allocator);
                defer pattern2.deinit();

                // Non-overlapping patterns are not valid neighbors
                if (!pattern1.overlaps(&pattern2, dx, dy)) {
                    continue;
                }

                if (pattern_to_idx.get(pattern_id2)) |pattern2_idx| {
                    neighbors.setValue(pattern2_idx, true);
                }
            }
            pattern_neigbors[i] = neighbors;
        }
        if (pattern_to_idx.get(pattern_id1)) |pattern1_idx| {
            propagator[pattern1_idx] = pattern_neigbors;
        }
    }
    return propagator;
}

pub fn propagate(self: *OverlapModel) !void {
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

            for (remove_idx.items) |remove_pattern_idx| {
                try self.remove_pattern(neig_idx, remove_pattern_idx);
            }
        }
    }
}

/// Single iteration of collapsing probabilities
///
/// Returns completion status
pub fn collapse(
    self: *OverlapModel,
    // Reusable allocated buffer to skip re-allocating each iteration
    random: std.Random,
) !BaseModel.Collapsed {
    var possible_pattern_weights = std.AutoArrayHashMap(usize, f32).init(self.allocator);
    try possible_pattern_weights.ensureTotalCapacity(self.pattern_weights.keys().len);

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
                const entropy = @as(f32, @log(@as(f32, @floatFromInt(weights_sum)))) - (weight_log_weight_sum / @as(f32, @floatFromInt(weights_sum)));

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
    while (argmin_possible_iter.next()) |pattern_idx| {
        const pattern_id = self.idx_to_pattern[pattern_idx];
        if (self.pattern_weights.get(pattern_id)) |weight| {
            possible_pattern_weights.putAssumeCapacity(pattern_idx, @as(f32, @floatFromInt(weight)));
        }
    }

    const random_idx = try random_pattern_idx(possible_pattern_weights, random.float(f32));

    var remove = try std.ArrayList(usize).initCapacity(self.allocator, possible_pattern_weights.keys().len);
    defer remove.deinit();
    var iter = self.base.possibilities[argmin_idx].iterator(.{});
    while (iter.next()) |pattern_idx| {
        if (pattern_idx != random_idx) {
            remove.appendAssumeCapacity(pattern_idx);
        }
    }

    std.debug.print("removing {} patterns from argmin_idx {}\n", .{ remove.items.len, argmin_idx });
    for (remove.items) |idx| {
        try self.remove_pattern(argmin_idx, idx);
    }
    return .failed;
}

pub fn collapse_ground(self: *OverlapModel) !void {
    if (self.ground_patterns.keys().len == 0) {
        return;
    }

    const ground_y: usize = self.base.output_cells_h - 1;
    var remove = try std.ArrayList(usize).initCapacity(
        self.allocator,
        self.pattern_weights.keys().len,
    );
    defer remove.deinit();

    for (0..self.base.output_cells_w) |x| {
        const ground_idx: usize = x + ground_y * self.base.output_cells_w;
        const ground_possibilities = self.base.possibilities[ground_idx];

        var ground_possiblities_iter = ground_possibilities.iterator(.{});
        while (ground_possiblities_iter.next()) |pattern_idx| {
            const pattern_id = self.idx_to_pattern[pattern_idx];
            if (self.ground_patterns.get(pattern_id) == null) {
                remove.appendAssumeCapacity(pattern_idx);
            }
        }

        for (remove.items) |remove_idx| {
            try self.remove_pattern(ground_idx, remove_idx);
        }
        remove.clearRetainingCapacity();
    }
}

/// Return the first index `x` where `sum of values [0 .. x - 1] > r`.
///
/// This is a form of doing weighted RNG, where the elements in `[0 .. x - 1]`
/// with higher value have a higher chance of being chosen, since greater
/// values have higher chance of being greater than `r`
fn random_pattern_idx(pattern_weights: *std.AutoArrayHashMap(usize, f32), r: f32) !usize {
    const sum = blk: {
        var s: f32 = 0.0;
        for (pattern_weights.values()) |v| {
            s += v;
        }
        if (s == 0) {
            @memset(pattern_weights.values(), 1.0);
            s = @as(f32, @floatFromInt(pattern_weights.values().len));
        }
        break :blk s;
    };

    var x: f32 = 0.0;
    var last_key: usize = undefined; // Keep track of the last valid key
    var count_iter = pattern_weights.iterator();
    while (count_iter.next()) |entry| {
        last_key = entry.key_ptr.*; // Store the key on every iteration
        const normalised = entry.value_ptr.* / sum;
        entry.value_ptr.* = normalised;
        x += normalised;
        if (r <= x) {
            return entry.key_ptr.*;
        }
    }
    return last_key;
}

pub fn remove_pattern(self: *OverlapModel, idx: usize, pattern_idx: usize) !void {
    // Already removed, break
    if (!self.base.possibilities[idx].isSet(pattern_idx)) {
        return;
    }

    self.base.possibilities[idx].unset(pattern_idx);
    self.base.possibilities_remaining[idx] -= 1;

    const pattern_id = self.idx_to_pattern[pattern_idx];
    if (self.pattern_weights.get(pattern_id)) |count| {
        self.base.weight_sums[idx] -= count;
    } else {
        std.log.err("remove_pattern: pattern_id: {d} has missing count\n", .{pattern_id});
    }

    if (self.weight_log_weights.get(pattern_id)) |weight_log_weight| {
        self.base.weight_log_weight_sums[idx] -= weight_log_weight;
    } else {
        std.log.err("remove_pattern: pattern_id: {d} has missing weight_log_weight\n", .{pattern_id});
    }

    const added_to_stack = self.base.on_stack.isSet(idx);
    if (!added_to_stack) {
        try self.base.propagation_stack.append(idx);
        self.base.on_stack.set(idx);
    }
}

pub fn clear(self: *OverlapModel) void {
    var weight_total: u32 = 0;
    var weight_log_weight_total: f32 = 0.0;
    var pattern_weights_iter = self.pattern_weights.iterator();
    while (pattern_weights_iter.next()) |entry| {
        const weight = entry.value_ptr.*;
        weight_total += weight;
        // Note: self.weight_log_weights is static and already calculated.
        // We can just sum its values.
        if (self.weight_log_weights.get(entry.key_ptr.*)) |wlw| {
            weight_log_weight_total += wlw;
        }
    }
    const num_patterns = self.pattern_weights.keys().len;
    self.base.clear(
        num_patterns,
        weight_total,
        weight_log_weight_total,
    );
}

/// Generates output pixels for the cell at `cell_idx`
pub fn render_cell(self: *OverlapModel, cell_idx: usize, output_pixel_idxs: []u32) !void {
    const cell_patterns = self.base.possibilities[cell_idx];
    const pattern_idx = cell_patterns.findFirstSet() orelse {
        return error.NoPossiblePatterns;
    };

    const pattern_id = self.idx_to_pattern[pattern_idx];
    const decoded = try self.pattern_encoder.decode(pattern_id, self.allocator);
    defer decoded.deinit();

    // NOTE: output pixel is the FIRST (top left) pixel of the pattern
    // since for overlap model cell_h and cell_w is 1x1
    const pixel_idx: u32 = decoded.values.items[0];
    output_pixel_idxs[cell_idx] = pixel_idx;
}

/// Output slice MUST be deinitialised
pub fn render_output(self: *OverlapModel) ![]u32 {
    const output_pixel_idxs: []u32 = try self.allocator.alloc(u32, self.base.output_cells_w * self.base.output_cells_h);
    for (0..self.base.output_cells_h) |y| {
        for (0..self.base.output_cells_w) |x| {
            const cell_idx = x + self.base.output_cells_w * y;
            try self.render_cell(cell_idx, output_pixel_idxs);
        }
    }
    return output_pixel_idxs;
}
