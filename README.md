# Codec2

## Usage

Create a Codec2 object (`Codec2::new(Codec2Mode::MODE_3200)`) then repeatedly obtain 8khz
16-bit raw audio samples and call `encode` to encode blocks of 160 samples into an 8-byte
(64-bit) compressed form, in order. On the receiving end, create a Codec2 object and repeatedly
call `decode` to decompress each chunk of 8 received bytes into the next 160 samples.

Complete example below.

Add this to your `Cargo.toml`:

```toml
[dependencies]
zerocopy = "0.3.0"
codec2 = "0"
```

main.rs:

```rust
use codec2::*;
use std::env::args;
use std::io::prelude::*;
use zerocopy::LayoutVerified;

fn main() -> std::io::Result<()> {
    if args().len() != 4 || (args().nth(1).unwrap() != "enc" && args().nth(1).unwrap() != "dec") {
        eprintln!("Usage: {} (enc|dec) inputfile outputfilename", args().nth(0).unwrap());
        eprintln!("Files should be raw 16-bit signed 8khz audio");
        return Ok(());
    }
    let mut fin = std::fs::File::open(args().nth(2).unwrap())?;
    let mut fout = std::fs::File::create(args().nth(3).unwrap())?;
    let mut c2 = Codec2::new(Codec2Mode::MODE_3200);
    let mut sampsbuf = vec![0; c2.samples_per_frame() * 2]; //u8 I/O buffer for 16-bit samples
    let mut packed = vec![0; (c2.bits_per_frame() + 7) / 8]; //u8 I/O buffer for encoded bits
    if args().nth(1).unwrap() == "enc" {
        while fin.read_exact(&mut sampsbuf).is_ok() {
            let samps = LayoutVerified::new_slice(&sampsbuf[..]).expect("bad alignment"); //assumes samples are in correct endianness
            c2.encode(&mut packed, &samps[..]);
            fout.write_all(&packed)?;
        }
    } else {
        while fin.read_exact(&mut packed).is_ok() {
            let mut samps = LayoutVerified::new_slice(&mut sampsbuf[..]).expect("bad alignment");
            c2.decode(&mut samps[..], &packed);
            fout.write_all(&sampsbuf)?;
        }
    }
    Ok(())
}
```
