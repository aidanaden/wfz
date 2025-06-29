const std = @import("std");
const Allocator = std.mem.Allocator;

const zigimg = @import("zigimg");
pub const Pixel = zigimg.color.Rgba32;

const Self = @This();

/// Image data
images: []zigimg.Image,
img_pixel_idxs: [][]u32,

/// Pixel to idx mappings
idx_to_pixel: std.ArrayList(Pixel),
pixel_to_idx: std.AutoArrayHashMap(Pixel, u32),
num_uniq_pixels: usize,

allocator: Allocator,

pub fn init(img_paths: [][]const u8, allocator: Allocator) !Self {
    var images = try allocator.alloc(zigimg.Image, img_paths.len);
    var img_pixel_idxs: [][]u32 = try allocator.alloc([]u32, img_paths.len);

    var pixel_to_idx = std.AutoArrayHashMap(Pixel, u32).init(allocator);
    var idx_to_pixel = std.ArrayList(Pixel).init(allocator);

    for (img_paths, 0..) |img_path, img_i| {
        var image = try zigimg.Image.fromFilePath(allocator, img_path);
        try image.convert(zigimg.PixelFormat.rgba32);

        const pixels = image.pixels.rgba32;
        var pixel_idxs: []u32 = try allocator.alloc(u32, pixels.len);

        var pixel_idx: u32 = 0;
        for (pixels, 0..) |pixel, pixel_i| {
            if (pixel_to_idx.get(pixel)) |idx| {
                pixel_idxs[pixel_i] = idx;
                continue;
            }
            // Handle newly found pixel value
            try pixel_to_idx.put(pixel, pixel_idx);
            try idx_to_pixel.append(pixel);
            pixel_idxs[pixel_i] = pixel_idx;
            pixel_idx += 1;
        }

        images[img_i] = image;
        img_pixel_idxs[img_i] = pixel_idxs;
    }

    return Self{
        .images = images,
        .img_pixel_idxs = img_pixel_idxs,
        .idx_to_pixel = idx_to_pixel,
        .pixel_to_idx = pixel_to_idx,
        .num_uniq_pixels = idx_to_pixel.items.len,
        .allocator = allocator,
    };
}

pub fn generate(self: *Self, output_w: usize, output_h: usize, output_pixel_idxs: []u32) !zigimg.Image {
    var output_pixels = try std.ArrayList(u8).initCapacity(self.allocator, output_w * output_h);
    defer output_pixels.deinit();

    // Build pixels from pixel idxs
    for (output_pixel_idxs) |pixel_idx| {
        const pixel = self.idx_to_pixel.items[pixel_idx];
        try output_pixels.appendSlice(std.mem.asBytes(&pixel));
    }

    // Generate image
    const pixel_fmt = self.images[0].pixelFormat();
    const output_img = try zigimg.Image.fromRawPixels(
        self.allocator,
        output_w,
        output_h,
        output_pixels.items[0..],
        pixel_fmt,
    );
    return output_img;
}

pub fn deinit(self: *Self) void {
    self.allocator.free(self.images);

    for (self.img_pixel_idxs) |pixel_idxs| {
        self.allocator.free(pixel_idxs);
    }
    self.idx_to_pixel.deinit();
    self.pixel_to_idx.deinit();
}
