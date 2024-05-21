use std::{
    error::Error,
    fs::File,
    io::{self, prelude::*, BufRead, BufReader},
    path::Path,
};

const MAGIC_MAX_LEN: usize = 6;
const GZ_MAGIC: [u8; 3] = [0x1f, 0x8b, 0x08];
const BZ_MAGIC: [u8; 3] = [0x42, 0x5a, 0x68];
const XZ_MAGIC: [u8; 6] = [0xfd, 0x37, 0x7a, 0x58, 0x5A, 0x00];
const BUFF_SIZE: usize = 512 * 1024;

fn magic_num<P: AsRef<Path> + Copy>(file_name: P) -> Result<[u8; MAGIC_MAX_LEN], Box<dyn Error>> {
    let mut buffer: [u8; MAGIC_MAX_LEN] = [0; MAGIC_MAX_LEN];
    let mut fp = File::open(file_name)?;
    let _ = fp.read(&mut buffer)?;
    Ok(buffer)
}

fn is_gzipped<P: AsRef<Path> + Copy>(file_name: P) -> Result<bool, Box<dyn Error>> {
    let buffer = magic_num(file_name)?;
    let gz_or_not =
        buffer[0] == GZ_MAGIC[0] && buffer[1] == GZ_MAGIC[1] && buffer[2] == GZ_MAGIC[2];
    Ok(gz_or_not
        || file_name
            .as_ref()
            .extension()
            .is_some_and(|ext| ext == "gz"))
}

fn is_bzipped<P: AsRef<Path> + Copy>(file_name: P) -> Result<bool, Box<dyn Error>> {
    let buffer = magic_num(file_name)?;
    let bz_or_not =
        buffer[0] == BZ_MAGIC[0] && buffer[1] == BZ_MAGIC[1] && buffer[2] == BZ_MAGIC[2];
    Ok(bz_or_not
        || file_name
            .as_ref()
            .extension()
            .is_some_and(|ext| ext == "bz2"))
}

fn is_xz<P: AsRef<Path> + Copy>(file_name: P) -> Result<bool, Box<dyn Error>> {
    let buffer = magic_num(file_name)?;
    let xz_or_not = buffer[0] == XZ_MAGIC[0]
        && buffer[1] == XZ_MAGIC[1]
        && buffer[2] == XZ_MAGIC[2]
        && buffer[3] == XZ_MAGIC[3]
        && buffer[4] == XZ_MAGIC[4]
        && buffer[5] == XZ_MAGIC[5];
    Ok(xz_or_not
        || file_name
            .as_ref()
            .extension()
            .is_some_and(|ext| ext == "xz"))
}

pub fn file_reader<P>(file_in: Option<P>) -> Result<Box<dyn BufRead + Send>, Box<dyn Error>>
where
    P: AsRef<Path> + Copy,
{
    if let Some(file_name) = file_in {
        let gz_flag = is_gzipped(file_name)?;
        let bz_flag = is_bzipped(file_name)?;
        let zx_flag = is_xz(file_name)?;

        let fp = File::open(file_name)?;

        if gz_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE,
                flate2::read::MultiGzDecoder::new(fp),
            )))
        } else if bz_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE,
                bzip2::read::MultiBzDecoder::new(fp),
            )))
        } else if zx_flag {
            Ok(Box::new(BufReader::with_capacity(
                BUFF_SIZE,
                xz2::read::XzDecoder::new_multi_decoder(fp),
            )))
        } else {
            Ok(Box::new(BufReader::with_capacity(BUFF_SIZE, fp)))
        }
    } else {
        if atty::is(atty::Stream::Stdin) {
            eprintln!("Error: stdin not detected");
            std::process::exit(1);
        }
        let fp = BufReader::new(io::stdin());
        Ok(Box::new(fp))
    }
}
