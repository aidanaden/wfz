const std = @import("std");
const math = std.math;
const Allocator = std.mem.Allocator;

/// Representation of an NxN tile of pixel values (typically pixel ids)
pub fn Pattern(V: type) type {
    return struct {
        N: usize,
        values: std.ArrayList(V),

        const Self = @This();

        pub fn init(start_x: usize, start_y: usize, N: usize, source: []V, source_h: usize, source_w: usize, allocator: Allocator) !Self {
            const size: usize = N * N;
            var values = try std.ArrayList(V).initCapacity(
                allocator,
                size,
            );
            for (0..N) |y_offset| {
                for (0..N) |x_offset| {
                    const source_x: usize = (start_x + x_offset) % source_w;
                    const source_y: usize = (start_y + y_offset) % source_h;
                    const source_idx = source_x + (source_y * source_w);
                    if (source_idx >= source.len) {
                        std.log.warn("invalid tile idx found: {}, skipping!\n", .{source_idx});
                        continue;
                    }
                    const value = source[source_idx];
                    const idx: usize = x_offset + y_offset * @as(usize, @intCast(N));
                    try values.insert(idx, value);
                }
            }
            return Self{
                .N = N,
                .values = values,
            };
        }

        pub fn deinit(self: *const Self) void {
            self.values.deinit();
        }

        pub const TransformMethod = enum {
            Rotate,
            Reflect,
        };

        pub fn transform(self: *const Self, method: TransformMethod) !Self {
            var values = try self.values.clone();
            values.clearRetainingCapacity();
            for (0..self.N) |og_y| {
                for (0..self.N) |og_x| {
                    const og_idx = og_x + og_y * self.N;
                    const x, const y = switch (method) {
                        .Reflect => .{ self.N - 1 - og_x, og_y },
                        .Rotate => .{ self.N - 1 - og_y, og_x },
                    };
                    const idx: usize = x + y * self.N;
                    try values.insert(og_idx, self.values.items[idx]);
                }
            }
            return Self{
                .N = self.N,
                .values = values,
            };
        }

        pub fn encode(self: *const Self, base: usize) BaseEncoder.EncodedId {
            return BaseEncoder.init(self.N, base).encode(self.values.items);
        }

        pub fn overlaps(self: *const Self, other: *const Self, dx: i32, dy: i32) bool {
            const i_N: i32 = @as(i32, @intCast(self.N));

            // Iteration 0: dx = dy = -2, xmin = ymin = 0, xmax = ymax = -2 + 3 = 1
            // Iteration 3: dx = dy = 1, xmin = ymin = 1, xmax = ymax = 3
            // Iteration 4: dx = dy = 2, xmin = ymin = 2, xmax = ymax = 3
            const xmin, const xmax = if (dx < 0) .{ 0, dx + i_N } else .{ dx, i_N };
            const ymin, const ymax = if (dy < 0) .{ 0, dy + i_N } else .{ dy, i_N };

            // If any value within the overlapping boundary
            // isnt the same, the patterns dont overlap
            var y = ymin;
            while (y < ymax) : (y += 1) {
                var x = xmin;
                while (x < xmax) : (x += 1) {
                    // Iteration 0: idx = xmin + ymin * 3 = 0 + 0 * 3 = 0
                    // Iteration 3.1: idx = xmin + ymin * 3 = 1 + 1 * 3 = 4
                    // Iteration 3.2: idx = xmin + ymin * 3 = 1 + 2 * 3 = 7
                    // Iteration 3.3: idx = xmin + ymin * 3 = 2 + 1 * 3 = 5
                    // Iteration 3.4: idx = xmin + ymin * 3 = 2 + 2 * 3 = 8
                    // Iteration 4: idx = xmin + ymin * 3 = 2 + 2 * 3 = 8
                    const idx: usize = @as(usize, @intCast(x + y * i_N));
                    // Calculate other pattern's matching `idx` value
                    // Iteration 0: other_idx = (0 - -2) + (0 - -2) * 3 = 2 + 2 * 3 = 8
                    // Iteration 3.1: other_idx = (1 - 1) + (1 - 1) * 3 = 0
                    // Iteration 3.3: other_idx = (1 - 1) + (2 - 1) * 3 = 3
                    // Iteration 3.2: other_idx = (2 - 1) + (1 - 1) * 3 = 1
                    // Iteration 3.4: other_idx = (2 - 1) + (2 - 1) * 3 = 4
                    // Iteration 4: other_idx = (2 - 2) + (2 - 2) * 3 = 0
                    const other_i_idx = (@as(i32, @intCast(x)) - dx) + ((@as(i32, @intCast(y)) - dy) * i_N);
                    const other_idx: usize = @as(usize, @intCast(other_i_idx));
                    if (self.values.items[idx] != other.values.items[other_idx]) {
                        return false;
                    }
                }
            }

            return true;
        }

        pub const BaseEncoder = struct {
            pub const EncodedId = usize;

            /// size = N x N
            size: usize,
            base: usize,

            pub fn init(N: usize, base: usize) BaseEncoder {
                const size: usize = N * N;
                return BaseEncoder{
                    .size = size,
                    .base = base,
                };
            }

            pub fn encode(self: *const BaseEncoder, values: []const V) EncodedId {
                var res: usize = 0;
                // Efficient calculation via horner's method
                for (values) |value| {
                    res = res * self.base + @as(usize, @intCast(value));
                }
                return res;
            }

            pub fn decode(self: *const BaseEncoder, id: EncodedId, allocator: Allocator) !Self {
                var values = try std.ArrayList(V).initCapacity(
                    allocator,
                    self.size,
                );
                try values.appendNTimes(0, self.size);

                var remainder = id;
                var i: usize = self.size - 1;
                while (remainder > 0) : ({
                    remainder = @divFloor(remainder, self.base);
                    if (i > 0) {
                        i -= 1;
                    }
                }) {
                    const digit = @mod(remainder, self.base);
                    values.items[i] = @as(u32, @intCast(digit));
                }
                return Self{
                    .N = @as(usize, @intFromFloat(@sqrt(@as(f32, @floatFromInt(self.size))))),
                    .values = values,
                };
            }
        };
    };
}

test "init, encode, decode" {
    const allocator = std.testing.allocator;

    const N = 3;
    const base = 256;
    const pattern = [_]u32{ 0, 0, 0, 0, 0, 0, 1, 2, 3 };
    const encoder = try Pattern(u32).BaseEncoder.init(N, base, allocator);

    // Encoding
    const encoded_id = encoder.encode(pattern[0..]);
    const expected_id = 66051;
    std.debug.print("encoded: {}, expected: {}\n", .{ encoded_id, expected_id });
    try std.testing.expect(encoded_id == expected_id);

    // Decoding
    const decoded_pattern = try encoder.decode(encoded_id, allocator);
    defer decoded_pattern.deinit();
    std.debug.print("decoded: {any}, pattern: {any}\n", .{ decoded_pattern.values, pattern });
    try std.testing.expect(std.mem.eql(u32, pattern[0..], decoded_pattern.values.items[0..]));
}
