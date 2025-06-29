const std = @import("std");
const Allocator = std.mem.Allocator;
const BaseModel = @import("../Base.zig");

const TileSymmetry = enum {
    /// 2 unique variants (2 rotations, 0 reflections)
    ///
    /// Rotations: 0, 90 degrees
    I,
    /// 1 unique variant (1 rotation, 0 reflections)
    ///
    /// Rotations: 0 degree
    X,
    /// 8 unique variants (4 rotations, 4 reflections)
    ///
    /// Rotations: 0, 90, 180, 270 degrees
    ///
    /// Reflections: 4 reflections
    L,
    /// 4 unique variants (4 rotations, 0 reflections)
    ///
    /// Rotations: 0, 90, 180, 270 degree
    T,
    /// 2 unique variants (2 reflections)
    ///
    /// Reflections: 2
    @"\"",
};

/// Represents an original non-permutated tile image
pub const SourceTile = struct {
    name: []const u8,
    symmetry: TileSymmetry,
    weight: f32 = 1.0,
    neighbors: []std.DynamicBitSet,
};

const TileSet = @This();

/// Original tiles data
source_tiles: []SourceTile,

/// Maps source tile name -> idx
source_idxs: std.StringArrayHashMap(u32),

allocator: Allocator,

pub fn fromJsonFile(filepath: []const u8, allocator: Allocator) !TileSet {
    std.debug.print("reading from json: {s}\n", .{filepath});
    const json_tileset = try readFromJsonFile(filepath, allocator);
    defer json_tileset.deinit();

    var source_tiles: []SourceTile = try allocator.alloc(SourceTile, json_tileset.value.tiles.len);

    var source_idxs = std.StringArrayHashMap(u32).init(allocator);
    try source_idxs.ensureTotalCapacity(json_tileset.value.tiles.len);
    for (json_tileset.value.tiles, 0..) |tile, idx| {
        source_idxs.putAssumeCapacity(tile.name, @as(u32, @intCast(idx)));
    }

    var neighbors = std.StringArrayHashMap([]std.DynamicBitSet).init(allocator);
    defer neighbors.deinit();

    // Iterate thru relations and build all possible
    // up/down/left/right neighbors for each tile
    for (json_tileset.value.neighbours) |relation| {
        var left_iter = std.mem.splitScalar(u8, relation.left, ' ');
        const left_name = left_iter.next() orelse {
            std.debug.print("invalid left relation found, missing left name, relation.left: {s}\n", .{relation.left});
            continue;
        };
        const left_idx = source_idxs.get(left_name) orelse {
            std.debug.print("invalid left tile found, relation.left: {s}\n", .{relation.left});
            continue;
        };

        // Set up neighbors for left tile
        var left_default_neighbors = try allocator.alloc(std.DynamicBitSet, BaseModel.Direction.len());
        for (0..BaseModel.Direction.len()) |i| {
            left_default_neighbors[i] = try std.DynamicBitSet.initEmpty(allocator, json_tileset.value.tiles.len);
        }
        const left_neighbors_entry = neighbors.getOrPutValue(left_name, left_default_neighbors) catch |err| {
            std.debug.print("failed to get left_neighbors, left_name: {s}, err: {any}\n", .{ left_name, err });
            continue;
        };
        const left_neighbors: []std.DynamicBitSet = left_neighbors_entry.value_ptr.*;

        // Extract direction of relation for left tile
        const left_dir = dir_blk: {
            const left_dir_idx_char = left_iter.next() orelse {
                // No idx found, default to right relation
                break :dir_blk BaseModel.Direction.Right;
            };
            // Invalid idx found, return null
            if (left_dir_idx_char.len > 1) {
                std.debug.print("invalid left_dir neighbor idx: {any}\n", .{left_dir_idx_char});
                break :dir_blk null;
            }
            switch (left_dir_idx_char[0]) {
                '1'...'3' => |raw_idx_char| {
                    const raw_idx: u8 = raw_idx_char - '0';
                    const raw_dir = @as(u2, @intCast(raw_idx));
                    const dir: BaseModel.Direction = @enumFromInt(raw_dir);
                    break :dir_blk dir;
                },
                else => |idx| {
                    std.debug.print("invalid left_dir neighbor direction idx found: {any}\n", .{idx});
                    break :dir_blk null;
                },
            }

            // All cases should be handled in switch statement
            comptime unreachable;
        };

        if (left_dir == null) {
            continue;
        }

        var right_iter = std.mem.splitScalar(u8, relation.right, ' ');
        const right_name = right_iter.next() orelse {
            std.debug.print("invalid relation found, relation.right: {s}\n", .{relation.right});
            continue;
        };
        const right_idx = source_idxs.get(right_name) orelse {
            std.debug.print("invalid right tile found, relation.right: {s}\n", .{relation.right});
            continue;
        };

        // Set up neighbors for right tile
        var right_default_neighbors = try allocator.alloc(std.DynamicBitSet, BaseModel.Direction.len());
        for (0..BaseModel.Direction.len()) |i| {
            right_default_neighbors[i] = try std.DynamicBitSet.initEmpty(allocator, json_tileset.value.tiles.len);
        }
        const right_neighbors_entry = neighbors.getOrPutValue(right_name, right_default_neighbors) catch |err| {
            std.debug.print("failed to get right_neighbors, right_name: {s}, err: {any}\n", .{ right_name, err });
            continue;
        };
        const right_neighbors: []std.DynamicBitSet = right_neighbors_entry.value_ptr.*;

        // Extract direction of relation for right tile
        const right_dir = dir_blk: {
            const right_dir_idx_char = left_iter.next() orelse {
                // No idx found, default to right relation
                break :dir_blk BaseModel.Direction.Right;
            };
            // Invalid idx found, return null
            if (right_dir_idx_char.len > 1) {
                std.debug.print("invalid right_dir neighbor idx: {any}\n", .{right_dir_idx_char});
                break :dir_blk null;
            }
            switch (right_dir_idx_char[0]) {
                '1'...'3' => |raw_idx_char| {
                    const raw_idx: u8 = raw_idx_char - '0';
                    // We inverse the direction idx since by default
                    // direction values are oriented from the pov of
                    // the left tile, so we inverse idx for right tile
                    const raw_idx_inverse: u8 = @mod(raw_idx + 2, 4);
                    const raw_dir = @as(u2, @intCast(raw_idx_inverse));
                    const dir: BaseModel.Direction = @enumFromInt(raw_dir);
                    break :dir_blk dir;
                },
                else => |idx| {
                    std.debug.print("invalid right_dir neighbor direction idx found: {any}\n", .{idx});
                    break :dir_blk null;
                },
            }

            // All cases should be handled in switch statement
            comptime unreachable;
        };

        if (right_dir == null) {
            continue;
        }

        // Add right tile idx as neighbor of left tile
        left_neighbors[@intFromEnum(left_dir.?)].setValue(@as(usize, @intCast(right_idx)), true);

        // Add left tile idx as neighbor of right tile
        right_neighbors[@intFromEnum(right_dir.?)].setValue(@as(usize, @intCast(left_idx)), true);
    }

    for (json_tileset.value.tiles, 0..) |json_tile, i| {
        const symmetry = std.meta.stringToEnum(TileSymmetry, json_tile.symmetry) orelse {
            continue;
        };
        const tile_neighbors: []std.DynamicBitSet = neighbors.get(json_tile.name) orelse {
            std.debug.print("tile has no neighbors, tile_name: {s}\n", .{json_tile.name});
            continue;
        };

        // std.debug.print("{s}: right neigbors: {d}, top neighbors: {d}, left neighbors: {d}, bottom neighbors: {d}\n", .{
        //     json_tile.name,
        //     tile_neighbors[0].count(),
        //     tile_neighbors[1].count(),
        //     tile_neighbors[2].count(),
        //     tile_neighbors[3].count(),
        // });

        for (0..BaseModel.Direction.len()) |dir_i| {
            const dir_neig = tile_neighbors[dir_i];
            if (dir_neig.count() == 0) {
                continue;
            }
            var dir_neig_names = try std.ArrayList([]const u8).initCapacity(allocator, dir_neig.count());
            defer dir_neig_names.deinit();
            var dir_neig_iter = dir_neig.iterator(.{});
            while (dir_neig_iter.next()) |idx| {
                const name = json_tileset.value.tiles[idx].name;
                dir_neig_names.appendAssumeCapacity(name);
            }
        }

        source_tiles[i] = SourceTile{
            .name = json_tile.name,
            .symmetry = symmetry,
            .weight = 1.0,
            .neighbors = tile_neighbors,
        };
    }
    return TileSet{
        .source_tiles = source_tiles,
        .source_idxs = source_idxs,
        .allocator = allocator,
    };
}

pub fn deinit(self: *TileSet) void {
    for (self.source_tiles) |*tile| {
        for (tile.neighbors) |*neig| {
            neig.deinit();
        }
        self.allocator.free(tile.neighbors);
    }
    self.allocator.free(self.source_tiles);
    self.source_idxs.deinit();
}

const JsonTileSet = struct {
    tiles: []JsonTile,
    neighbours: []JsonTileRelation,
};

const JsonTile = struct {
    name: []const u8,
    symmetry: []const u8,
};

const JsonTileRelation = struct {
    left: []const u8,
    right: []const u8,
};

fn readFromJsonFile(filepath: []const u8, allocator: Allocator) !std.json.Parsed(JsonTileSet) {
    const json_file = try std.fs.cwd().openFile(filepath, .{});
    defer json_file.close();
    var buf_reader = std.io.bufferedReader(json_file.reader());
    var reader = buf_reader.reader();

    const stat = try json_file.stat();
    const lines = try reader.readAllAlloc(allocator, stat.size);
    const json = try std.json.parseFromSlice(JsonTileSet, allocator, lines, .{});
    return json;
}

pub fn print_tile_neighbors(self: *const TileSet) !void {
    for (self.source_tiles) |source_tile| {
        const tile_neighbors = source_tile.neighbors;
        std.debug.print("{s} neighbors:\n", .{source_tile.name});
        for (0..BaseModel.Direction.len()) |dir_i| {
            const dir_neig = tile_neighbors[dir_i];
            if (dir_neig.count() == 0) {
                continue;
            }
            var dir_neig_names = try std.ArrayList([]const u8).initCapacity(self.allocator, dir_neig.count());
            defer dir_neig_names.deinit();
            var dir_neig_iter = dir_neig.iterator(.{});
            while (dir_neig_iter.next()) |idx| {
                const name = self.source_tiles[idx].name;
                dir_neig_names.appendAssumeCapacity(name);
            }

            const dir_neig_log = try std.mem.join(self.allocator, ", ", dir_neig_names.items);
            defer self.allocator.free(dir_neig_log);
            std.debug.print("{s} ({d}): {s}\n", .{ @tagName(@as(BaseModel.Direction, @enumFromInt(dir_i))), dir_neig_names.items.len, dir_neig_log });
        }
        std.debug.print("\n", .{});
    }
}

// pub fn fromXmlFile(filename: []const u8) !TileSet {}

test "read from json file" {
    const test_fname = "./src/assets/tilesets/castle.json";
    const json_tileset = try readFromJsonFile(test_fname, std.testing.allocator);
    defer json_tileset.deinit();
    std.debug.print("num tiles: {}, num neighbors: {}\n", .{ json_tileset.value.tiles.len, json_tileset.value.neighbours.len });
}
