pub const tracy = @import("tracy");
pub const tracy_impl = @import("tracy_impl");
pub const tracy_options: tracy.Options = .{
    .default_callstack_depth = 8,
    .verbose = true,
};

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const Collapsed = enum {
    success,
    contradiction,
    failed,
};

/// --- SHARED DATA ---
/// Width of output in pixels
output_w: usize,

/// Width of output in cells
output_cells_w: usize,

/// Height of output in pixels
output_h: usize,

/// Height of output in cells
output_cells_h: usize,

/// Whether output should be repeating
periodic: bool,

/// The "Wave Function": A bitset of possible choices for each **CELL**,
/// updated whenever `propagate` is called
possibilities: []std.DynamicBitSet,

/// Cached data for optimised entropy calculation
possibilities_remaining: []usize,
weight_sums: []f32,
weight_log_weight_sums: []f32,

/// Propagation tracking
/// Used to track cells that must be updated
propagation_stack: std.ArrayList(usize),

/// Track if cell is on stack to be propagated
/// Used for optimization purposes
on_stack: std.DynamicBitSet,

allocator: Allocator,

const BaseModel = @This();

pub const Direction = enum(u2) {
    // 0 degree
    Right,
    // 90 degree
    Up,
    // 180 degree
    Left,
    // 270 degree
    Down,
    pub fn point(self: Direction) struct { i32, i32 } {
        return switch (self) {
            Direction.Right => .{ 1, 0 },
            Direction.Up => .{ 0, 1 },
            Direction.Left => .{ -1, 0 },
            Direction.Down => .{ 0, -1 },
        };
    }
    pub fn len() usize {
        return std.meta.tags(Direction).len;
    }
};

pub const InitParams = struct {
    N: usize,
    periodic: bool,
    output_w: usize,
    output_h: usize,
    symmetry: usize,
    ground: bool,
    name: []const u8,
    img_path: [:0]const u8,
    output_fname: []const u8,
};

pub fn init(
    params: *const InitParams,
    // Width of output image in terms of number of cells
    output_cells_w: usize,
    // Height of output image in terms of number of cells
    output_cells_h: usize,
    num_possibilities: usize,
    weight_total: f32,
    weight_log_weight_total: f32,
    allocator: Allocator,
) !BaseModel {
    var possibilities: []std.DynamicBitSet = try allocator.alloc(std.DynamicBitSet, output_cells_h * output_cells_w);
    for (0..output_cells_h) |y| {
        for (0..output_cells_w) |x| {
            const neighbors = try std.DynamicBitSet.initFull(allocator, num_possibilities);
            const idx = x + y * output_cells_w;
            possibilities[idx] = neighbors;
        }
    }

    // Store weight (pattern count) of all possible patterns for each cell
    const weight_sums = try allocator.alloc(f32, output_cells_h * output_cells_w);

    // Store weight * log weight of all possible patterns for each cell
    const weight_log_weight_sums = try allocator.alloc(f32, output_cells_h * output_cells_w);

    // Store number of possible patterns remaining per cell
    const possibilities_remaining = try allocator.alloc(usize, output_cells_h * output_cells_w);

    @memset(weight_sums, weight_total);
    @memset(weight_log_weight_sums, weight_log_weight_total);
    @memset(possibilities_remaining, num_possibilities);

    const propagation_stack = std.ArrayList(usize).init(allocator);
    const on_stack = try std.DynamicBitSet.initEmpty(allocator, output_cells_h * output_cells_w);

    return BaseModel{
        .output_w = params.output_w,
        .output_cells_w = output_cells_w,
        .output_h = params.output_h,
        .output_cells_h = output_cells_h,
        .periodic = params.periodic,
        .possibilities = possibilities,
        .weight_sums = weight_sums,
        .weight_log_weight_sums = weight_log_weight_sums,
        .possibilities_remaining = possibilities_remaining,
        .propagation_stack = propagation_stack,
        .on_stack = on_stack,
        .allocator = allocator,
    };
}

pub fn deinit(self: *BaseModel) void {

    // Possibilities
    for (self.possibilities) |*neighbors| {
        neighbors.deinit();
    }
    self.allocator.free(self.possibilities);

    // Cached weights/entropy values
    self.allocator.free(self.weight_sums);
    self.allocator.free(self.weight_log_weight_sums);
    self.allocator.free(self.possibilities_remaining);

    // Propagation stack
    self.propagation_stack.deinit();
    self.on_stack.deinit();
}

pub fn clear(self: *BaseModel, num_patterns: usize, weight_total: f32, weight_log_weight_total: f32) void {
    for (0..self.possibilities.len) |i| {
        self.possibilities[i].setRangeValue(.{ .start = 0, .end = num_patterns }, true);
    }

    @memset(self.weight_sums, weight_total);
    // Initialise all weight_log_weight_sums to total weight_log_weight
    @memset(self.weight_log_weight_sums, weight_log_weight_total);
    // Initialise patterns remaining to total number of patterns  (all patterns possible at the start)
    @memset(self.possibilities_remaining, num_patterns);

    self.propagation_stack.clearRetainingCapacity();
    self.on_stack.setRangeValue(
        .{
            .start = 0,
            .end = self.on_stack.capacity() - 1,
        },
        false,
    );
}
