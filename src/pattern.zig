const std = @import("std");
const Allocator = std.mem.Allocator;

pub fn Pattern(V: type) type {
    return struct {
        N: usize,
        values: std.ArrayList(V),
        allocator: Allocator,

        const Self = @This();

        pub fn init(start_x: usize, start_y: usize, N: usize, source: []V, source_h: usize, source_w: usize, allocator: Allocator) !Self {
            const size: usize = N * N;
            var values = try std.ArrayList(V).initCapacity(
                allocator,
                @as(usize, @intCast(size)),
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
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *const Self) void {
            self.values.deinit();
        }

        pub fn rotate(self: *const Self) !Self {
            var rotated_values = try std.ArrayList(usize).initCapacity(
                self.allocator,
                self.values.len,
            );
            for (0..self.N) |y| {
                for (0..self.N) |x| {
                    const idx = x + y * self.N;
                    const rotated_x = self.N - 1 - y;
                    const rotated_y = x;
                    const rotated_idx = rotated_x + rotated_y * self.N;
                    try rotated_values.insert(idx, self.values[rotated_idx]);
                }
            }
            return Self{
                .N = self.N,
                .values = rotated_values,
                .allocator = self.allocator,
            };
        }

        pub const BaseEncoder = struct {
            const EncodedId = usize;
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
                    res = res * self.base + @as(usize, @intCast(@as(u32, @bitCast(value))));
                }
                return res;
            }

            const DecodeError = error{
                InvalidSize,
            };
            pub fn decode(self: *const BaseEncoder, id: EncodedId, decoded: []usize) DecodeError!void {
                // Invalid size, decoded slice provided is not NxN
                if (decoded.len != self.size) {
                    return DecodeError.InvalidSize;
                }
                var remainder = id;
                var i: usize = self.size - 1;
                while (remainder > 0) : ({
                    remainder = @divFloor(remainder, self.base);
                    i -= 1;
                }) {
                    const digit = @mod(remainder, self.base);
                    decoded[i] = digit;
                }
            }
        };
    };
}

test "InitEncodeDecode" {
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
    var decoded: [9]u32 = std.mem.zeroes([9]u32);
    try encoder.decode(encoded_id, decoded[0..]);
    std.debug.print("decoded: {any}, pattern: {any}\n", .{ decoded, pattern });
    try std.testing.expect(std.mem.eql(u32, pattern[0..], decoded[0..]));
}
