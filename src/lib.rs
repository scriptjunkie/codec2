//! Codec 2 is an open source (LGPL 2.1) low bit rate speech codec.
//!
//! This is a zero dependencies pure rust port of the original [http://rowetel.com/codec2.html](http://rowetel.com/codec2.html)
//! Currently 3200 and 2400 bitrates encoding and decoding are implemented.
//!
//! # Basic Usage
//!
//! Create a Codec2 object e.g. `Codec2::new(Codec2Mode::MODE_3200)` then repeatedly obtain raw 8khz
//! 16-bit audio samples and call `encode` to encode blocks of `samples_per_frame()` samples into
//! `bits_per_frame()` compressed bits, in order. On the receiving end, create a Codec2 object and
//! repeatedly call `decode` to decompress each chunk of received bytes into the next samples.
//!
//! Complete example below. This example uses zerocopy to interpret the `&[i16]` slices as `&[u8]`
//! for I/O.
//!
//! Cargo.toml:
//!
//! ```toml
//! [dependencies]
//! zerocopy = "0.3.0"
//! codec2 = "0"
//! ```
//!
//! main.rs:
//!
//! ```rust
//! use codec2::*;
//! use std::env::args;
//! use std::io::prelude::*;
//! use zerocopy::AsBytes;
//!
//! fn main() -> std::io::Result<()> {
//!     if args().len() != 4 || (args().nth(1).unwrap() != "enc" && args().nth(1).unwrap() != "dec") {
//!         eprintln!("Usage: {} (enc|dec) inputfile outputfilename", args().nth(0).unwrap());
//!         eprintln!("Files should be raw 16-bit signed 8khz audio");
//!         return Ok(());
//!     }
//!     let mut fin = std::fs::File::open(args().nth(2).unwrap())?;
//!     let mut fout = std::fs::File::create(args().nth(3).unwrap())?;
//!     let mut c2 = Codec2::new(Codec2Mode::MODE_3200);
//!     let mut samps = vec![0; c2.samples_per_frame()]; //u16 I/O buffer
//!     let mut packed = vec![0; (c2.bits_per_frame() + 7) / 8]; //u8 I/O buffer for encoded bits
//!     if args().nth(1).unwrap() == "enc" {
//!         while fin.read_exact(samps.as_bytes_mut()).is_ok() {
//!             c2.encode(&mut packed, &samps[..]);
//!             fout.write_all(&packed)?;
//!         }
//!     } else {
//!         while fin.read_exact(&mut packed).is_ok() {
//!             c2.decode(&mut samps[..], &packed);
//!             fout.write_all(samps.as_bytes())?;
//!         }
//!     }
//!     Ok(())
//! }
//! ```

#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![warn(trivial_numeric_casts)]

#[cfg(test)]
mod tests {
    #[test]
    fn encode_test() {
        let mut c2 = crate::Codec2::new(crate::Codec2Mode::MODE_3200);
        //160 samples
        let spf = c2.samples_per_frame();
        let mut samps = vec![0; spf * 2];
        let bpf = c2.bits_per_frame();
        println!(
            "c2 is n_samp {:?} m_pitch {:?} Sn.len {:?} bpf {} spf {}",
            c2.internal.n_samp,
            c2.internal.m_pitch,
            c2.internal.Sn.len(),
            bpf,
            spf
        );
        let mut outbuf = vec![0; (bpf + 7) / 8];
        c2.encode(&mut outbuf, &samps);
        assert_eq!(&outbuf[..], &[0xC0, 0, 0x6A, 0x43, 0x9C, 0xE4, 0x21, 8][..]);
        let outbuforig = outbuf.clone();
        println!("encoded {:X?}", outbuf);
        c2.encode(&mut outbuf, &samps);
        assert_eq!(&outbuf[..], &[0x81, 0, 9, 0x43, 0x9C, 0xE4, 0x21, 8][..]);
        println!("encoded {:X?}", outbuf);
        c2.encode(&mut outbuf, &samps);
        assert_eq!(&outbuf[..], &[1, 0, 9, 0x43, 0x9C, 0xE4, 0x21, 8][..]);
        println!("encoded {:X?}", outbuf);
        c2.encode(&mut outbuf, &samps);
        assert_eq!(&outbuf[..], &[1, 0, 9, 0x43, 0x9C, 0xE4, 0x21, 8][..]);
        println!("encoded {:X?}", outbuf);

        c2.decode(&mut samps, &outbuforig);
        println!("decoded {:X?}", samps);
        for samp in &samps {
            assert!((*samp).abs() < 0x100);
        }
    }
}

mod kiss_fft;
mod nlp;
use nlp::FDMDV_OS_TAPS_16K;
mod quantise;
use crate::quantise::*;
mod codebook;
use crate::codebook::*;
mod codebookd;
use crate::codebookd::*;
mod codebookge;
use crate::codebookge::*;
const WO_BITS: i32 = 7;
const WO_E_BITS: u32 = 8;
const LSPD_SCALAR_INDEXES: usize = 10;
const LSP_SCALAR_INDEXES: usize = 10;
use std::f64::consts::PI;

const N_S: f32 = 0.01; //  internal proc frame length in secs
const TW_S: f32 = 0.005; //  trapezoidal synth window overlap
const MAX_AMP: usize = 160; //  maximum number of harmonics
const TWO_PI: f32 = 6.283185307; //  mathematical constant
const FFT_ENC: usize = 512; //  size of FFT used for encoder
const FFT_DEC: usize = 512; //  size of FFT used in decoder
const V_THRESH: f32 = 6.0; //  voicing threshold in dB
const LPC_ORD: usize = 10; //  LPC order
                           //  Pitch estimation constants
const M_PITCH_S: f32 = 0.0400; //  pitch analysis window in s
const P_MIN_S: f32 = 0.0025; //  minimum pitch period in s
const P_MAX_S: f32 = 0.0200; //  maximum pitch period in s
mod inner {
    use crate::*;
    #[derive(Clone, Debug)]
    pub struct C2const {
        pub Fs: i32,       //  sample rate of this instance
        pub n_samp: usize, //  number of samples per 10ms frame at Fs
        //        pub max_amp: i32,   //  maximum number of harmonics
        pub m_pitch: usize, //  pitch estimation window size in samples
        pub p_min: i32,     //  minimum pitch period in samples
        pub p_max: i32,     //  maximum pitch period in samples
        pub Wo_min: f32,
        pub Wo_max: f32,
        pub nw: usize, //  analysis window size in samples
        pub tw: usize, //  trapezoidal synthesis window overlap
    }
    impl C2const {
        pub fn new(Fs: i32, framelength_s: f32) -> Self {
            Self {
                Fs: Fs,
                n_samp: ((Fs as f32) * framelength_s).round() as usize,
                //                max_amp: ((Fs as f32) * P_MAX_S / 2.0).floor() as i32,
                p_min: ((Fs as f32) * P_MIN_S).floor() as i32,
                p_max: ((Fs as f32) * P_MAX_S).floor() as i32,
                m_pitch: ((Fs as f32) * M_PITCH_S).floor() as usize,
                Wo_min: TWO_PI / ((Fs as f32) * P_MAX_S).floor(),
                Wo_max: TWO_PI / ((Fs as f32) * P_MIN_S).floor(),
                nw: 279,
                tw: ((Fs as f32) * TW_S) as usize,
            }
        }

        /*---------------------------------------------------------------------------*\

          FUNCTION....: dft_speech
          AUTHOR......: David Rowe, conversion by Matt Weeks
          DATE CREATED: 27/5/94

          Finds the DFT of the current speech input speech frame.

        \*---------------------------------------------------------------------------*/
        // TODO: we can either go for a faster FFT using fftr and some stack usage
        // or we can reduce stack usage to almost zero on STM32 by switching to fft_inplace
        pub fn dft_speech(
            &self,
            fft_fwd_cfg: &codec2_fft_cfg,
            Sw: &mut [COMP],
            Sn: &[f32],
            w: &[f32],
        ) {
            let m_pitch = self.m_pitch;
            let nw = self.nw;

            for i in 0..FFT_ENC {
                Sw[i].r = 0.0;
                Sw[i].i = 0.0;
            }

            /* Centre analysis window on time axis, we need to arrange input
            to FFT this way to make FFT phases correct */

            //  move 2nd half to start of FFT input vector

            for i in 0..nw / 2 {
                Sw[i].r = Sn[i + m_pitch / 2] * w[i + m_pitch / 2];
            }

            //  move 1st half to end of FFT input vector

            for i in 0..nw / 2 {
                Sw[FFT_ENC - nw / 2 + i].r =
                    Sn[i + m_pitch / 2 - nw / 2] * w[i + m_pitch / 2 - nw / 2];
            }

            codec2_fft_inplace(fft_fwd_cfg, Sw);
        }
    }
    //  describes each codebook
    #[derive(Clone, Debug)]
    pub struct lsp_codebook {
        pub k: usize,           //  dimension of vector
        pub log2m: i32,         //  number of bits in m
        pub cb: &'static [f32], //  The elements
    }
    #[derive(Copy, Clone)]
    pub struct kiss_fft_cpx {
        pub r: kiss_fft_scalar,
        pub i: kiss_fft_scalar,
    }
    impl kiss_fft_cpx {
        pub fn new() -> Self {
            Self { r: 0.0, i: 0.0 }
        }
        pub fn kf_cexp(phase: f64) -> Self {
            Self {
                r: phase.cos() as f32,
                i: phase.sin() as f32,
            }
        }
    }
    impl std::fmt::Debug for kiss_fft_cpx {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_fmt(format_args!("{},{}", self.r, self.i))
        }
    }
    pub type COMP = kiss_fft_cpx;
    #[derive(Clone, Debug)]
    pub struct kiss_fft_state {
        pub nfft: usize,
        pub inverse: i32,
        pub factors: [usize; 2 * MAXFACTORS],
        pub twiddles: Vec<kiss_fft_cpx>,
    }
    impl kiss_fft_state {
        pub fn new(nfft: usize, inverse_fft: i32) -> Self {
            let mut twiddles = Vec::with_capacity(nfft);
            let mut factors = [0; 2 * MAXFACTORS];
            for i in 0..nfft {
                let mut phase = -2.0 * PI * (i as f64) / (nfft as f64);
                if inverse_fft != 0 {
                    phase *= -1.0;
                }
                twiddles.push(kiss_fft_cpx::kf_cexp(phase));
            }
            let mut n = nfft;
            let mut p = 4;
            let floor_sqrt = (n as f64).sqrt().floor() as usize;
            let mut idx = 0;
            // factor out powers of 4, powers of 2, then any remaining primes
            loop {
                while (n % p) != 0 {
                    match p {
                        4 => p = 2,
                        2 => p = 3,
                        _ => p += 2,
                    }
                    if p > floor_sqrt {
                        p = n; //  no more factors, skip to end
                    }
                }
                n /= p;
                factors[idx] = p;
                idx += 1;
                factors[idx] = n;
                idx += 1;
                if n <= 1 {
                    break;
                } //do..while (n > 1)
            }
            Self {
                nfft: nfft,
                inverse: inverse_fft,
                factors: factors,
                twiddles: twiddles,
            }
        }
    }
    #[derive(Clone, Debug)]
    pub struct kiss_fftr_state {
        pub substate: kiss_fft_cfg,
        pub tmpbuf: Vec<kiss_fft_cpx>,
        pub super_twiddles: Vec<kiss_fft_cpx>,
    }
    impl kiss_fftr_state {
        pub fn new(mut nfft: usize, inverse_fft: i32) -> Self {
            nfft = nfft / 2;
            let mut res = Self {
                substate: kiss_fft_state::new(nfft, inverse_fft),
                tmpbuf: vec![kiss_fft_cpx::new(); nfft],
                super_twiddles: vec![kiss_fft_cpx::new(); nfft / 2],
            };
            for i in 0..nfft / 2 {
                let mut phase =
                    -3.14159265358979323846264338327 * ((i as f64 + 1.0) / nfft as f64 + 0.5);
                if inverse_fft != 0 {
                    phase *= -1.0;
                }
                res.super_twiddles[i] = kiss_fft_cpx::kf_cexp(phase);
            }
            res
        }
    }
    #[derive(Clone, Debug)]
    pub struct NLP {
        pub Fs: i32, //  sample rate in Hz
        pub m: usize,
        pub w: [f32; PMAX_M / DEC], //  DFT window
        pub sq: [f32; PMAX_M],      //  squared speech samples
        pub mem_x: f32,
        pub mem_y: f32,               //  memory for notch filter
        pub mem_fir: [f32; NLP_NTAP], //  decimation FIR filter memory
        pub fft_cfg: codec2_fft_cfg,  //  kiss FFT config
        pub Sn16k: Vec<f32>,          //  Fs=16kHz input speech vector
                                      //    FILE         *f,
    }
    impl NLP {
        pub fn new(c2const: &C2const) -> Self {
            let (m, vsize) = if c2const.Fs == 16000 {
                (
                    c2const.m_pitch / 2,
                    FDMDV_OS_TAPS_16K as usize + c2const.n_samp,
                )
            } else {
                (c2const.m_pitch, 0)
            };
            let mut w = [0.0; PMAX_M / DEC];
            for i in 0..m / DEC {
                w[i] =
                    0.5 - 0.5 * (2.0 * PI as f32 * i as f32 / (m as f32 / DEC as f32 - 1.0)).cos();
            }
            Self {
                Fs: c2const.Fs, //  sample rate in Hz
                m: m,
                w: w,              //  DFT window
                sq: [0.0; PMAX_M], //  squared speech samples
                mem_x: 0.0,
                mem_y: 0.0,                                   //  memory for notch filter
                mem_fir: [0.0; NLP_NTAP],                     //  decimation FIR filter memory
                fft_cfg: kiss_fft_state::new(PE_FFT_SIZE, 0), //  kiss FFT config
                Sn16k: vec![0.0; vsize],                      //  Fs=16kHz input speech vector
            }
        }
    }
    #[derive(Clone, Debug, Copy)]
    pub struct MODEL {
        pub Wo: f32,                 //  fundamental frequency estimate in radians
        pub L: usize,                //  number of harmonics
        pub A: [f32; MAX_AMP + 1],   //  amplitiude of each harmonic
        pub phi: [f32; MAX_AMP + 1], //  phase of each harmonic
        pub voiced: i32,             //  non-zero if this frame is voiced
    }
    impl MODEL {
        pub fn new(p_max: f32) -> Self {
            let wo = TWO_PI / p_max;
            Self {
                Wo: wo,                               //  fundamental frequency estimate in radians
                L: (PI / wo as f64).floor() as usize, //  number of harmonics
                A: [0.0; MAX_AMP + 1],                //  amplitiude of each harmonic
                phi: [0.0; MAX_AMP + 1],              //  phase of each harmonic
                voiced: 0,                            //  non-zero if this frame is voiced
            }
        }
    }
}
use inner::*;
type kiss_fft_scalar = f32;
type kiss_fft_cfg = kiss_fft_state;

/* e.g. an fft of length 128 has 4 factors
as far as kissfft is concerned
4*4*4*2
*/
const PE_FFT_SIZE: usize = 512;
const MAXFACTORS: usize = 32;

type codec2_fft_cfg = kiss_fft_state;
type codec2_fftr_cfg = kiss_fftr_state;
type kiss_fftr_cfg = kiss_fftr_state;
fn codec2_fft(cfg: &kiss_fft_cfg, fin: &[kiss_fft_cpx], fout: &mut [kiss_fft_cpx]) {
    kiss_fft::kiss_fft(cfg, fin, fout);
}
fn codec2_fftr(cfg: &mut codec2_fftr_cfg, inp: &[f32], out: &mut [kiss_fft_cpx]) {
    kiss_fft::kiss_fftr(cfg, inp, out);
}
fn codec2_fftri(cfg: &mut codec2_fftr_cfg, inp: &[kiss_fft_cpx], out: &mut [f32]) {
    kiss_fft::kiss_fftri(cfg, inp, out);
}

fn codec2_fft_inplace(cfg: &codec2_fft_cfg, inout: &mut [kiss_fft_cpx]) {
    let mut in_ = [kiss_fft_cpx::new(); 512];
    // decide whether to use the local stack based buffer for in
    // or to allow kiss_fft to allocate RAM
    // second part is just to play safe since first method
    // is much faster and uses less RAM
    if cfg.nfft <= 512 {
        in_[..cfg.nfft].copy_from_slice(&inout[..cfg.nfft]);
        kiss_fft::kiss_fft(cfg, &in_, inout);
    } else {
        kiss_fft::kiss_fft(cfg, &inout.to_vec(), inout);
    }
}

const PMAX_M: usize = 320;
const DEC: usize = 5;
const NLP_NTAP: usize = 48;
const BPF_N: usize = 101;
const LPCPF_BETA: f32 = 0.2;
const LPCPF_GAMMA: f32 = 0.5;

#[derive(Clone, Debug)]
struct Codec2Internal {
    mode: Codec2Mode,
    c2const: C2const,
    Fs: i32,
    n_samp: usize,
    m_pitch: usize,
    fft_fwd_cfg: codec2_fft_cfg,   //  forward FFT config
    fftr_fwd_cfg: codec2_fftr_cfg, //  forward real FFT config
    w: Vec<f32>,                   //  [m_pitch] time domain hamming window
    W: [f32; FFT_ENC],             //  DFT of w[]
    Pn: Vec<f32>,                  //  [2*n_samp] trapezoidal synthesis window
    bpf_buf: Vec<f32>,             //  buffer for band pass filter
    Sn: Vec<f32>,                  //  [m_pitch] input speech
    hpf_states: [f32; 2],          //  high pass filter states
    nlp: NLP,                      //  pitch predictor states
    gray: i32,                     //  non-zero for gray encoding

    fftr_inv_cfg: codec2_fftr_cfg, //  inverse FFT config
    Sn_: Vec<f32>,                 //  [2*n_samp] synthesised output speech
    ex_phase: f32,                 //  excitation model phase track
    bg_est: f32,                   //  background noise estimate for post filter
    prev_f0_enc: f32,              //  previous frame's f0    estimate
    prev_model_dec: MODEL,         //  previous frame's model parameters
    prev_lsps_dec: [f32; LPC_ORD], //  previous frame's LSPs
    prev_e_dec: f32,               //  previous frame's LPC energy

    lpc_pf: i32,     //  LPC post filter on
    bass_boost: i32, //  LPC post filter bass boost
    beta: f32,       //  LPC post filter parameters
    gamma: f32,

    xq_enc: [f32; 2], //  joint pitch and energy VQ states
    xq_dec: [f32; 2],

    smoothing: i32, //  enable smoothing for channels with errors
    se: f32,        //  running sum of squared error
    nse: u32,       //  number of terms in sum
    post_filter_en: i32,
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: make_synthesis_window
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 11/5/94

  Init function that generates the trapezoidal (Parzen) sythesis window.

\*---------------------------------------------------------------------------*/
fn make_synthesis_window(c2const: &C2const, Pn: &mut [f32]) {
    let n_samp = c2const.n_samp;
    let tw = c2const.tw;

    //  Generate Parzen window in time domain
    for i in 0..n_samp / 2 - tw {
        Pn[i] = 0.0;
    }
    let mut win = 0.0;
    for i in n_samp / 2 - tw..n_samp / 2 + tw {
        Pn[i] = win;
        win += 1.0 / (2.0 * tw as f32);
    }
    for i in n_samp / 2 + tw..3 * n_samp / 2 - tw {
        Pn[i] = 1.0;
    }
    win = 1.0;
    for i in 3 * n_samp / 2 - tw..3 * n_samp / 2 + tw {
        Pn[i] = win;
        win -= 1.0 / (2.0 * tw as f32);
    }
    for i in 3 * n_samp / 2 + tw..2 * n_samp {
        Pn[i] = 0.0;
    }
}

fn make_analysis_window(
    c2const: &C2const,
    fft_fwd_cfg: &codec2_fft_cfg,
    w: &mut [f32],
    W: &mut [f32],
) {
    let mut wshift = [kiss_fft_cpx::new(); FFT_ENC];
    let m_pitch = c2const.m_pitch;
    let nw = c2const.nw;

    /*
       Generate Hamming window centered on M-sample pitch analysis window

    0            M/2           M-1
    |-------------|-------------|
          |-------|-------|
              nw samples

       All our analysis/synthsis is centred on the M/2 sample.
    */

    let mut m = 0.0;
    for i in 0..m_pitch / 2 - nw / 2 {
        w[i] = 0.0;
    }
    let mut j = 0;
    for i in m_pitch / 2 - nw / 2..m_pitch / 2 + nw / 2 {
        w[i] = 0.5 - 0.5 * (TWO_PI * (j as f32) / ((nw as f32) - 1.0)).cos();
        m += w[i] * w[i];
        j += 1;
    }
    for i in m_pitch / 2 + nw / 2..m_pitch {
        w[i] = 0.0;
    }

    /* Normalise - makes freq domain amplitude estimation straight
    forward */

    m = 1.0 / (m * FFT_ENC as f32).sqrt();
    for i in 0..m_pitch {
        w[i] *= m;
    }

    /*
       Generate DFT of analysis window, used for later processing.  Note
       we modulo FFT_ENC shift the time domain window w[], this makes the
       imaginary part of the DFT W[] equal to zero as the shifted w[] is
       even about the n=0 time axis if nw is odd.  Having the imag part
       of the DFT W[] makes computation easier.

       0                      FFT_ENC-1
       |-------------------------|

        ----\               /----
             \             /
              \           /          <- shifted version of window w[n]
               \         /
                \       /
                 -------

       |---------|     |---------|
         nw/2              nw/2
    */

    let mut temp = [kiss_fft_cpx::new(); FFT_ENC];

    for i in 0..FFT_ENC {
        wshift[i].r = 0.0;
        wshift[i].i = 0.0;
    }
    for i in 0..(nw / 2) {
        wshift[i].r = w[i + (m_pitch) / 2];
    }
    let mut j = m_pitch / 2 - nw / 2;
    for i in FFT_ENC - nw / 2..FFT_ENC {
        wshift[i].r = w[j];
        j += 1;
    }

    codec2_fft(fft_fwd_cfg, &wshift[..], &mut temp);

    /*
        Re-arrange W[] to be symmetrical about FFT_ENC/2.  Makes later
        analysis convenient.

     Before:


       0                 FFT_ENC-1
       |----------|---------|
       __                   _
         \                 /
          \_______________/

     After:

       0                 FFT_ENC-1
       |----------|---------|
                 ___
                /   \
       ________/     \_______

    */

    for i in 0..FFT_ENC / 2 {
        W[i] = temp[i + FFT_ENC / 2].r;
        W[i + FFT_ENC / 2] = temp[i].r;
    }
}

const WordSize: usize = 8;
const IndexMask: usize = 0x7;
const ShiftRight: usize = 3;
fn pack(
    bitArray: &mut [u8],  //  The output bit string.
    bitIndex: &mut usize, //  Index into the string in BITS, not bytes.
    field: i32,           //  The bit field to be packed.
    fieldWidth: u32,      //  Width of the field in BITS, not bytes.
) {
    pack_natural_or_gray(bitArray, bitIndex, field, fieldWidth, 1);
}

/** Pack a bit field into a bit string, encoding the field in Gray code.
 *
 * The output is an array of unsigned char data. The fields are efficiently
 * packed into the bit string. The Gray coding is a naive attempt to reduce
 * the effect of single-bit errors, we expect to do a better job as the
 * codec develops.
 *
 * This code would be simpler if it just set one bit at a time in the string,
 * but would hit the same cache line more often. I'm not sure the complexity
 * gains us anything here.
 *
 * Although field is currently of int type rather than unsigned for
 * compatibility with the rest of the code, indices are always expected to
 * be >= 0.
 */
fn pack_natural_or_gray(
    bitArray: &mut [u8],  //  The output bit string.
    bitIndex: &mut usize, //  Index into the string in BITS, not bytes.
    mut field: i32,       //  The bit field to be packed.
    mut fieldWidth: u32,  //  Width of the field in BITS, not bytes.
    gray: u32,            //  non-zero for gray coding
) {
    if gray > 0 {
        //  Convert the field to Gray code
        field = (field >> 1) ^ field;
    }

    while fieldWidth != 0 {
        let bI = *bitIndex;
        let bitsLeft = (WordSize - (bI & IndexMask)) as u32;
        let sliceWidth = if bitsLeft < fieldWidth {
            bitsLeft
        } else {
            fieldWidth
        };
        let wordIndex = bI >> ShiftRight;

        bitArray[wordIndex] |=
            (((field >> (fieldWidth - sliceWidth)) << (bitsLeft - sliceWidth)) & 0xFF) as u8;

        *bitIndex = bI + sliceWidth as usize;
        fieldWidth -= sliceWidth;
    }
}

/** Unpack a field from a bit string, converting from Gray code to binary.
 *
 */
fn unpack(
    bitArray: &[u8],    //  The input bit string.
    bitIndex: &mut u32, //  Index into the string in BITS, not bytes.
    fieldWidth: u32,    //  Width of the field in BITS, not bytes.
) -> i32 {
    unpack_natural_or_gray(bitArray, bitIndex, fieldWidth, 1)
}

/** Unpack a field from a bit string, to binary, optionally using
 * natural or Gray code.
 *
 */
fn unpack_natural_or_gray(
    bitArray: &[u8],     //  The input bit string.
    bitIndex: &mut u32,  //  Index into the string in BITS, not bytes.
    mut fieldWidth: u32, //  Width of the field in BITS, not bytes.
    gray: u32,           //  non-zero for Gray coding
) -> i32 {
    let mut field = 0;
    while fieldWidth != 0 {
        let bI = *bitIndex;
        let bitsLeft = WordSize as u32 - (bI & IndexMask as u32);
        let sliceWidth = if bitsLeft < fieldWidth {
            bitsLeft
        } else {
            fieldWidth
        };

        field |= ((bitArray[(bI >> ShiftRight) as usize] as usize >> (bitsLeft - sliceWidth))
            & ((1 << sliceWidth) - 1))
            << (fieldWidth - sliceWidth);

        *bitIndex = bI + sliceWidth;
        fieldWidth -= sliceWidth;
    }

    if gray != 0 {
        //  Convert from Gray code to binary. Works for maximum 8-bit fields.
        let mut t = field ^ (field >> 8);
        t ^= t >> 4;
        t ^= t >> 2;
        t ^= t >> 1;
        t as i32
    } else {
        field as i32
    }
}

/// Codec mode (bitrate). Currently only MODE_3200 and MODE_2400 are implemented
#[derive(Clone, Debug, Copy)]
pub enum Codec2Mode {
    MODE_3200,
    MODE_2400,
    //MODE_1600,
    //MODE_1400,
    //MODE_1300,
    //MODE_1200,
    //MODE_700C,
    //MODE_450,
    //MODE_450PWB,
}
use Codec2Mode::*;

/// A Codec2 object for encoding or decoding audio
#[derive(Clone, Debug)]
pub struct Codec2 {
    internal: Codec2Internal,
}

impl Codec2 {
    /// Creates a new Codec2 object suitable for encoding or decoding audio
    pub fn new(mode: Codec2Mode) -> Self {
        let c2const = C2const::new(8000, N_S);
        let n_samp = c2const.n_samp;
        let m_pitch = c2const.m_pitch;
        let mut c2 = Self {
            internal: Codec2Internal {
                mode: mode,
                Fs: c2const.Fs,
                n_samp: n_samp,
                m_pitch: m_pitch,
                fft_fwd_cfg: kiss_fft_state::new(FFT_ENC, 0), //  forward FFT config
                fftr_fwd_cfg: kiss_fftr_state::new(FFT_ENC, 0), //  forward real FFT config
                w: vec![0.0; m_pitch], //  [m_pitch] time domain hamming window
                W: [0.0; FFT_ENC],     //  DFT of w[]
                Pn: vec![0.0; 2 * n_samp], //  [2*n_samp] trapezoidal synthesis window
                bpf_buf: vec![0.0; BPF_N + 4 * n_samp], //  buffer for band pass filter
                Sn: vec![1.0; m_pitch], //  [m_pitch] input speech
                hpf_states: [0.0; 2],  //  high pass filter states

                nlp: NLP::new(&c2const), //  pitch predictor states
                gray: 1,                 //  non-zero for gray encoding

                fftr_inv_cfg: kiss_fftr_state::new(FFT_DEC, 1), //  inverse FFT config
                Sn_: vec![0.0; m_pitch], //  [2*n_samp] synthesised output speech
                prev_f0_enc: 1.0 / P_MAX_S,
                bg_est: 0.0,
                ex_phase: 0.0,
                prev_model_dec: MODEL::new(c2const.p_max as f32),
                c2const: c2const,
                prev_lsps_dec: [0.0; LPC_ORD],
                prev_e_dec: 1.0,
                lpc_pf: 1,        //  LPC post filter on
                bass_boost: 1,    //  LPC post filter bass boost
                beta: LPCPF_BETA, //  LPC post filter parameters
                gamma: LPCPF_GAMMA,
                xq_enc: [0.0; 2], //  joint pitch and energy VQ states
                xq_dec: [0.0; 2],

                smoothing: 0, //  enable smoothing for channels with errors
                se: 0.0,
                nse: 0,
                post_filter_en: 1,
            },
        };
        for i in 0..LPC_ORD {
            c2.internal.prev_lsps_dec[i] = (i as f64 * PI / (LPC_ORD as f64 + 1.0)) as f32;
        }
        make_analysis_window(
            &c2.internal.c2const,
            &c2.internal.fft_fwd_cfg,
            &mut c2.internal.w,
            &mut c2.internal.W,
        );
        make_synthesis_window(&c2.internal.c2const, &mut c2.internal.Pn);
        c2
    }

    /// The number of bits in an encoded (compressed) frame; 64 for the 3200 bitrate, 48 for 2400.
    pub fn bits_per_frame(&self) -> usize {
        if let MODE_3200 = self.internal.mode {
            64
        } else {
            48
        }
    }

    /// How many samples an encoded (compressed) frame represents; generally 160 (20ms of speech).
    pub fn samples_per_frame(&self) -> usize {
        160
    }

    /// Encodes speech samples at current bitrate into `bits_per_frame()`/8 rounded up output bytes.
    /// For MODE_3200, this is 64 bits or 8 bytes, for MODE_2400, it's 48 bits (6 bytes).
    pub fn encode(&mut self, bits: &mut [u8], speech: &[i16]) {
        if let MODE_3200 = self.internal.mode {
            self.codec2_encode_3200(bits, speech)
        } else {
            self.codec2_encode_2400(bits, speech)
        }
    }

    /// Decodes `bits_per_frame()` compressed bits into `samples_per_frame()` speech samples.
    /// For MODE_3200, the input is 64 bits or 8 bytes, for MODE_2400, it's 48 bits (6 bytes).
    pub fn decode(&mut self, speech: &mut [i16], bits: &[u8]) {
        if let MODE_3200 = self.internal.mode {
            self.codec2_decode_3200(speech, bits)
        } else {
            self.codec2_decode_2400(speech, bits)
        }
    }

    /*---------------------------------------------------------------------------*\

      FUNCTION....: codec2_encode_3200
      AUTHOR......: David Rowe, conversion by Matt Weeks
      DATE CREATED: 13 Sep 2012

      Encodes 160 speech samples (20ms of speech) into 64 bits.

      The codec2 algorithm actually operates internally on 10ms (80
      sample) frames, so we run the encoding algorithm twice.  On the
      first frame we just send the voicing bits.  On the second frame we
      send all model parameters.  Compared to 2400 we use a larger number
      of bits for the LSPs and non-VQ pitch and energy.

      The bit allocation is:

        Parameter                      bits/frame
        --------------------------------------
        Harmonic magnitudes (LSPs)     50
        Pitch (Wo)                      7
        Energy                          5
        Voicing (10ms update)           2
        TOTAL                          64

    \*---------------------------------------------------------------------------*/
    /// Encodes 160 speech samples (20ms of speech) into 64 bits.
    fn codec2_encode_3200(&mut self, bits: &mut [u8], speech: &[i16]) {
        let mut model = MODEL::new(self.internal.c2const.p_max as f32);
        let mut ak = [0.0; LPC_ORD + 1]; //f32
        let mut lsps = [0.0; LPC_ORD]; //f32
        let mut lspd_indexes = [0; LPC_ORD];
        let mut nbit = 0;

        let nbyte = (self.bits_per_frame() + 7) / 8;
        for i in 0..nbyte {
            bits[i] = 0;
        }

        //  first 10ms analysis frame - we just want voicing
        self.analyse_one_frame(&mut model, speech);
        pack(bits, &mut nbit, model.voiced, 1);

        //  second 10ms analysis frame
        self.analyse_one_frame(&mut model, &speech[self.internal.n_samp..]);
        pack(bits, &mut nbit, model.voiced, 1);
        let Wo_index = encode_Wo(&self.internal.c2const, model.Wo, WO_BITS);
        pack(bits, &mut nbit, Wo_index, WO_BITS as u32); //1+1+7 = 9 bits

        let e = speech_to_uq_lsps(
            &mut lsps,
            &mut ak,
            &self.internal.Sn,
            &self.internal.w,
            self.internal.m_pitch,
            LPC_ORD,
        );
        let e_index = encode_energy(e, E_BITS);
        pack(bits, &mut nbit, e_index, E_BITS as u32); //9+5 = 14 bits

        encode_lspds_scalar(&mut lspd_indexes, &lsps, LPC_ORD);
        for i in 0..LSPD_SCALAR_INDEXES {
            pack(bits, &mut nbit, lspd_indexes[i], lspd_bits(i) as u32);
        }
    }

    /*---------------------------------------------------------------------------*\

      FUNCTION....: codec2_decode_3200
      AUTHOR......: David Rowe, conversion by Matt Weeks
      DATE CREATED: 13 Sep 2012

      Decodes a frame of 64 bits into 160 samples (20ms) of speech.

    \*---------------------------------------------------------------------------*/
    /// Decodes a frame of 64 bits into 160 samples (20ms) of speech.
    fn codec2_decode_3200(&mut self, speech: &mut [i16], bits: &[u8]) {
        let mut ak = [[0.0; LPC_ORD + 1]; 2];
        let mut nbit = 0;
        let mut Aw = [COMP::new(); FFT_ENC];

        let mut model = [MODEL::new(self.internal.c2const.p_max as f32); 2];

        //  unpack bits from channel ------------------------------------

        /* this will partially fill the model params for the 2 x 10ms
        frames */

        model[0].voiced = unpack(bits, &mut nbit, 1);
        model[1].voiced = unpack(bits, &mut nbit, 1);

        let Wo_index = unpack(bits, &mut nbit, WO_BITS as u32);
        model[1].Wo = decode_Wo(&self.internal.c2const, Wo_index, WO_BITS);
        model[1].L = (PI / model[1].Wo as f64) as usize;

        let mut e = [0.0; 2];
        let e_index = unpack(bits, &mut nbit, E_BITS as u32);
        e[1] = decode_energy(e_index, E_BITS);

        let mut lspd_indexes = [0; LPC_ORD];
        for i in 0..LSPD_SCALAR_INDEXES {
            lspd_indexes[i] = unpack(bits, &mut nbit, lspd_bits(i) as u32) as usize;
        }
        let mut lsps = [[0.0; LPC_ORD]; 2];
        decode_lspds_scalar(&mut lsps[1][0..], &lspd_indexes, LPC_ORD);

        //  interpolate ------------------------------------------------

        /* Wo and energy are sampled every 20ms, so we interpolate just 1
        10ms frame between 20ms samples */

        let (model0, model1) = model.split_at_mut(1);
        interp_Wo(
            &mut model0[0],
            &self.internal.prev_model_dec,
            &model1[0],
            self.internal.c2const.Wo_min,
        );
        e[0] = interp_energy(self.internal.prev_e_dec, e[1]);

        /* LSPs are sampled every 20ms so we interpolate the frame in
        between, then recover spectral amplitudes */

        let (lsps0, lsps1) = lsps.split_at_mut(1);
        interpolate_lsp_ver2(
            &mut lsps0[0][0..],
            &self.internal.prev_lsps_dec,
            &mut lsps1[0][0..],
            0.5,
            LPC_ORD,
        );

        for i in 0..2 {
            lsp_to_lpc(&lsps[i][0..], &mut ak[i][0..], LPC_ORD);
            let mut snr = 0.0;
            aks_to_M2(
                &mut self.internal.fftr_fwd_cfg,
                &ak[i][..],
                LPC_ORD,
                &mut model[i],
                e[i],
                &mut snr,
                0,
                0,
                self.internal.lpc_pf,
                self.internal.bass_boost,
                self.internal.beta,
                self.internal.gamma,
                &mut Aw,
            );
            apply_lpc_correction(&mut model[i]);
            self.synthesise_one_frame(
                &mut speech[self.internal.n_samp * i..],
                &mut model[i],
                &mut Aw,
                1.0,
            );
        }

        //  update memories for next frame ----------------------------

        self.internal.prev_model_dec = model[1];
        self.internal.prev_e_dec = e[1];
        for i in 0..LPC_ORD {
            self.internal.prev_lsps_dec[i] = lsps[1][i];
        }
    }

    /*---------------------------------------------------------------------------*\

      FUNCTION....: codec2_encode_2400
      AUTHOR......: David Rowe
      DATE CREATED: 21/8/2010

      Encodes 160 speech samples (20ms of speech) into 48 bits.

      The codec2 algorithm actually operates internally on 10ms (80
      sample) frames, so we run the encoding algorithm twice.  On the
      first frame we just send the voicing bit.  On the second frame we
      send all model parameters.

      The bit allocation is:

        Parameter                      bits/frame
        --------------------------------------
        Harmonic magnitudes (LSPs)     36
        Joint VQ of Energy and Wo       8
        Voicing (10ms update)           2
        Spare                           2
        TOTAL                          48

    \*---------------------------------------------------------------------------*/
    /// Encodes 160 speech samples (20ms of speech) into 48 bits.
    fn codec2_encode_2400(&mut self, bits: &mut [u8], speech: &[i16]) {
        let mut model = MODEL::new(self.internal.c2const.p_max as f32);
        let mut ak = [0.0; LPC_ORD + 1]; //f32
        let mut lsps = [0.0; LPC_ORD]; //f32
        let mut lsp_indexes = [0; LPC_ORD];
        let mut nbit = 0;

        for i in 0..(self.bits_per_frame() + 7) / 8 {
            bits[i] = 0;
        }

        //  first 10ms analysis frame - we just want voicing
        self.analyse_one_frame(&mut model, speech);
        pack(bits, &mut nbit, model.voiced, 1);

        //  second 10ms analysis frame
        self.analyse_one_frame(&mut model, &speech[self.internal.n_samp..]);
        pack(bits, &mut nbit, model.voiced, 1);
        let e = speech_to_uq_lsps(
            &mut lsps,
            &mut ak,
            &self.internal.Sn,
            &self.internal.w,
            self.internal.m_pitch,
            LPC_ORD,
        );

        let WoE_index = encode_WoE(&model, e, &mut self.internal.xq_enc);
        pack(bits, &mut nbit, WoE_index, WO_E_BITS);

        encode_lsps_scalar(&mut lsp_indexes, &lsps, LPC_ORD);
        for i in 0..LSP_SCALAR_INDEXES {
            pack(bits, &mut nbit, lsp_indexes[i], lsp_bits(i));
        }
        let spare = 0;
        pack(bits, &mut nbit, spare, 2);

        //assert(nbit == (unsigned)codec2_bits_per_frame(c2));
    }

    /*---------------------------------------------------------------------------*\

      FUNCTION....: codec2_decode_2400
      AUTHOR......: David Rowe
      DATE CREATED: 21/8/2010

      Decodes frames of 48 bits into 160 samples (20ms) of speech.

    \*---------------------------------------------------------------------------*/
    fn codec2_decode_2400(&mut self, speech: &mut [i16], bits: &[u8]) {
        let mut model = [MODEL::new(self.internal.c2const.p_max as f32); 2];
        let mut lsp_indexes = [0; LPC_ORD];
        let mut lsps = [[0.0; LPC_ORD]; 2];
        let mut e = [0.0; 2];
        let mut snr = 0.0;
        let mut ak = [[0.0; LPC_ORD + 1]; 2];
        let mut nbit = 0;
        let mut Aw = [COMP::new(); FFT_ENC];

        //assert(c2 != NULL);

        /* unpack bits from channel ------------------------------------*/

        /* this will partially fill the model params for the 2 x 10ms
        frames */
        model[0].voiced = unpack(bits, &mut nbit, 1);
        model[1].voiced = unpack(bits, &mut nbit, 1);

        let WoE_index = unpack(bits, &mut nbit, WO_E_BITS) as usize;
        decode_WoE(
            &self.internal.c2const,
            &mut model[1],
            &mut e[1],
            &mut self.internal.xq_dec,
            WoE_index,
        );

        for i in 0..LSP_SCALAR_INDEXES {
            lsp_indexes[i] = unpack(bits, &mut nbit, lsp_bits(i)) as usize;
        }
        decode_lsps_scalar(&mut lsps[1][0..], &lsp_indexes, LPC_ORD);
        check_lsp_order(&mut lsps[1][0..], LPC_ORD);
        bw_expand_lsps(&mut lsps[1][0..], LPC_ORD, 50.0, 100.0);

        /* interpolate ------------------------------------------------*/

        /* Wo and energy are sampled every 20ms, so we interpolate just 1
        10ms frame between 20ms samples */

        let (model0, model1) = model.split_at_mut(1);
        interp_Wo(
            &mut model0[0],
            &self.internal.prev_model_dec,
            &model1[0],
            self.internal.c2const.Wo_min,
        );
        e[0] = interp_energy(self.internal.prev_e_dec, e[1]);

        /* LSPs are sampled every 20ms so we interpolate the frame in
        between, then recover spectral amplitudes */

        let (lsps0, lsps1) = lsps.split_at_mut(1);
        interpolate_lsp_ver2(
            &mut lsps0[0][0..],
            &self.internal.prev_lsps_dec,
            &mut lsps1[0][0..],
            0.5,
            LPC_ORD,
        );

        for i in 0..2 {
            lsp_to_lpc(&lsps[i][0..], &mut ak[i][0..], LPC_ORD);
            aks_to_M2(
                &mut self.internal.fftr_fwd_cfg,
                &ak[i][..],
                LPC_ORD,
                &mut model[i],
                e[i],
                &mut snr,
                0,
                0,
                self.internal.lpc_pf,
                self.internal.bass_boost,
                self.internal.beta,
                self.internal.gamma,
                &mut Aw,
            );
            apply_lpc_correction(&mut model[i]);
            self.synthesise_one_frame(
                &mut speech[self.internal.n_samp * i..],
                &mut model[i],
                &mut Aw,
                1.0,
            );
        }

        /* update memories for next frame ----------------------------*/
        self.internal.prev_model_dec = model[1];
        self.internal.prev_e_dec = e[1];
        for i in 0..LPC_ORD {
            self.internal.prev_lsps_dec[i] = lsps[1][i];
        }
    }

    /*---------------------------------------------------------------------------* \

      FUNCTION....: analyse_one_frame()
      AUTHOR......: David Rowe, conversion by Matt Weeks
      DATE CREATED: 23/8/2010

      Extract sinusoidal model parameters from 80 speech samples (10ms of
      speech).

    \*---------------------------------------------------------------------------*/
    fn analyse_one_frame(&mut self, model: &mut MODEL, speech: &[i16]) {
        let mut Sw = [COMP::new(); FFT_ENC];
        let n_samp = self.internal.n_samp;
        let m_pitch = self.internal.m_pitch;

        //  Read input speech
        for i in 0..m_pitch - n_samp {
            self.internal.Sn[i] = self.internal.Sn[i + n_samp];
        }
        for i in 0..n_samp {
            self.internal.Sn[i + m_pitch - n_samp] = speech[i].into();
        }
        self.internal.c2const.dft_speech(
            &self.internal.fft_fwd_cfg,
            &mut Sw,
            &self.internal.Sn,
            &self.internal.w,
        );

        //  Estimate pitch
        let mut pitch = 0.0;
        nlp::nlp(
            &mut self.internal.nlp,
            &self.internal.Sn,
            n_samp,
            &mut pitch,
            &Sw,
            &self.internal.W,
            &mut self.internal.prev_f0_enc,
        );
        model.Wo = TWO_PI / pitch;
        model.L = (PI / model.Wo as f64) as usize;

        //  estimate model parameters
        two_stage_pitch_refinement(&self.internal.c2const, model, &Sw);

        //  estimate phases when doing ML experiments
        estimate_amplitudes(model, &Sw, &self.internal.W, 0);
        est_voicing_mbe(&self.internal.c2const, model, &Sw, &self.internal.W);
    }

    /*---------------------------------------------------------------------------* \

      FUNCTION....: synthesise_one_frame()
      AUTHOR......: David Rowe, conversion by Matt Weeks
      DATE CREATED: 23/8/2010

      Synthesise 80 speech samples (10ms) from model parameters.

    \*---------------------------------------------------------------------------*/
    fn synthesise_one_frame(
        &mut self,
        speech: &mut [i16],
        model: &mut MODEL,
        Aw: &[COMP],
        gain: f32,
    ) {
        //  LPC based phase synthesis
        let mut H = [COMP::new(); MAX_AMP + 1];
        sample_phase(model, &mut H, Aw);
        phase_synth_zero_order(
            self.internal.n_samp,
            model,
            &mut self.internal.ex_phase,
            &mut H,
        );

        postfilter(model, &mut self.internal.bg_est);
        synthesise(
            self.internal.n_samp,
            &mut self.internal.fftr_inv_cfg,
            &mut self.internal.Sn_,
            model,
            &self.internal.Pn,
            true,
        );

        for i in 0..self.internal.n_samp {
            self.internal.Sn_[i] *= gain;
        }

        ear_protection(&mut self.internal.Sn_, self.internal.n_samp);

        for i in 0..self.internal.n_samp {
            if self.internal.Sn_[i] > 32767.0 {
                speech[i] = 32767;
            } else if self.internal.Sn_[i] < -32767.0 {
                speech[i] = -32767;
            } else {
                speech[i] = self.internal.Sn_[i] as i16;
            }
        }
    }
}

/*---------------------------------------------------------------------------* \

  FUNCTION....: ear_protection()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: Nov 7 2012

  Limits output level to protect ears when there are bit errors or the input
  is overdriven.  This doesn't correct or mask bit errors, just reduces the
  worst of their damage.

\*---------------------------------------------------------------------------*/
fn ear_protection(in_out: &mut [f32], n: usize) {
    //  find maximum sample in frame

    let mut max_sample = 0.0;
    for i in 0..n {
        if in_out[i] > max_sample {
            max_sample = in_out[i];
        }
    }

    //  determine how far above set point

    let over = max_sample / 30000.0;

    /* If we are x dB over set point we reduce level by 2x dB, this
    attenuates major excursions in amplitude (likely to be caused
    by bit errors) more than smaller ones */

    if over > 1.0 {
        let gain = 1.0 / (over * over);
        for i in 0..n {
            in_out[i] *= gain;
        }
    }
}

/*---------------------------------------------------------------------------*\

  est_voicing_mbe()

  Returns the error of the MBE cost function for a fiven F0.

  Note: I think a lot of the operations below can be simplified as
  W[].i = 0 and has been normalised such that den always equals 1.

\*---------------------------------------------------------------------------*/
fn est_voicing_mbe(c2const: &C2const, model: &mut MODEL, Sw: &[COMP], W: &[f32]) -> f32 {
    let mut Am = COMP::new(); // amplitude sample for this band
    let mut Ew = COMP::new();

    let l_1000hz = (model.L as f32 * 1000.0 / ((c2const.Fs / 2) as f32)) as usize;
    let mut sig = 1E-4;
    for l in 1..l_1000hz + 1 {
        sig += model.A[l] * model.A[l];
    }

    let Wo = model.Wo;
    let mut error = 1E-4; //accumulated error between original and synthesised

    //  Just test across the harmonics in the first 1000 Hz

    for l in 1..l_1000hz + 1 {
        Am.r = 0.0;
        Am.i = 0.0;
        let mut den = 0.0; //denominator of Am expression
        let al = ((l as f32 - 0.5) * Wo * (FFT_ENC as f32) / TWO_PI).ceil() as usize;
        let bl = ((l as f32 + 0.5) * Wo * (FFT_ENC as f32) / TWO_PI).ceil() as usize;

        //  Estimate amplitude of harmonic assuming harmonic is totally voiced

        // centers Hw[] about current harmonic
        let offset =
            (FFT_ENC as f32 / 2.0 - (l as f32) * Wo * (FFT_ENC as f32) / TWO_PI + 0.5) as usize;
        for m in al..bl {
            Am.r += Sw[m].r * W[offset + m];
            Am.i += Sw[m].i * W[offset + m];
            den += W[offset + m] * W[offset + m];
        }

        Am.r = Am.r / den;
        Am.i = Am.i / den;

        //  Determine error between estimated harmonic and original

        for m in al..bl {
            Ew.r = Sw[m].r - Am.r * W[offset + m];
            Ew.i = Sw[m].i - Am.i * W[offset + m];
            error += Ew.r * Ew.r;
            error += Ew.i * Ew.i;
        }
    }

    let snr = 10.0 * (sig / error).log10();
    if snr > V_THRESH {
        model.voiced = 1;
    } else {
        model.voiced = 0;
    }
    //  post processing, helps clean up some voicing errors ------------------

    /*
       Determine the ratio of low freqency to high frequency energy,
       voiced speech tends to be dominated by low frequency energy,
       unvoiced by high frequency. This measure can be used to
       determine if we have made any gross errors.
    */

    let l_2000hz = (model.L as f32 * 2000.0 / (c2const.Fs as f32 / 2.0)) as usize;
    let l_4000hz = (model.L as f32 * 4000.0 / (c2const.Fs as f32 / 2.0)) as usize;
    let mut ehigh = 1E-4;
    let mut elow = ehigh;
    for l in 1..l_2000hz + 1 {
        elow += model.A[l] * model.A[l];
    }
    for l in l_2000hz..l_4000hz + 1 {
        ehigh += model.A[l] * model.A[l];
    }
    let eratio = 10.0 * (elow / ehigh).log10();

    /* Look for Type 1 errors, strongly V speech that has been
    accidentally declared UV */

    if model.voiced == 0 {
        if eratio > 10.0 {
            model.voiced = 1;
        }
    }
    /* Look for Type 2 errors, strongly UV speech that has been
    accidentally declared V */

    if model.voiced == 1 {
        if eratio < -10.0 {
            model.voiced = 0;
        }
        /* A common source of Type 2 errors is the pitch estimator
        gives a low (50Hz) estimate for UV speech, which gives a
        good match with noise due to the close harmoonic spacing.
        These errors are much more common than people with 50Hz3
        pitch, so we have just a small eratio threshold. */

        let sixty = 60.0 * TWO_PI / (c2const.Fs as f32);
        if (eratio < -4.0) && (model.Wo <= sixty) {
            model.voiced = 0;
        }
    }

    return snr;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: estimate_amplitudes
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 27/5/94

  Estimates the complex amplitudes of the harmonics.

\*---------------------------------------------------------------------------*/
fn estimate_amplitudes(model: &mut MODEL, Sw: &[COMP], _W: &[f32], est_phase: i32) {
    let r = TWO_PI / (FFT_ENC as f32);
    let one_on_r = 1.0 / r;

    for m in 1..model.L + 1 {
        //  Estimate ampltude of harmonic

        let mut den = 0.0; // denominator of amplitude expression
                           // bounds of current harmonic
        let am = ((m as f32 - 0.5) * model.Wo * one_on_r + 0.5) as usize;
        let bm = ((m as f32 + 0.5) * model.Wo * one_on_r + 0.5) as usize;

        for i in am..bm {
            den += Sw[i].r * Sw[i].r + Sw[i].i * Sw[i].i;
        }

        model.A[m] = den.sqrt();

        if est_phase != 0 {
            let b = (m as f32 * model.Wo / r + 0.5) as usize; //  DFT bin of centre of current harmonic

            /* Estimate phase of harmonic, this is expensive in CPU for
            embedded devicesso we make it an option */

            model.phi[m] = Sw[b].i.atan2(Sw[b].r);
        }
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: two_stage_pitch_refinement
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 27/5/94

  Refines the current pitch estimate using the harmonic sum pitch
  estimation technique.

\*---------------------------------------------------------------------------*/
fn two_stage_pitch_refinement(c2const: &C2const, model: &mut MODEL, Sw: &[COMP]) {
    //  Coarse refinement
    //  pitch refinment minimum, maximum and step
    let mut pmax = TWO_PI / model.Wo + 5.0;
    let mut pmin = TWO_PI / model.Wo - 5.0;
    let mut pstep = 1.0;
    hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

    //  Fine refinement

    pmax = TWO_PI / model.Wo + 1.0;
    pmin = TWO_PI / model.Wo - 1.0;
    pstep = 0.25;
    hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

    //  Limit range

    if model.Wo < TWO_PI / (c2const.p_max as f32) {
        model.Wo = TWO_PI / (c2const.p_max as f32);
    }
    if model.Wo > TWO_PI / (c2const.p_min as f32) {
        model.Wo = TWO_PI / (c2const.p_min as f32);
    }

    model.L = (PI / model.Wo as f64).floor() as usize;

    //  trap occasional round off issues with floorf()
    if model.Wo * model.L as f32 >= 0.95 * PI as f32 {
        model.L -= 1;
    }
    //  assert(model.Wo*model.L < PI);
}

/*---------------------------------------------------------------------------*\

 FUNCTION....: hs_pitch_refinement
 AUTHOR......: David Rowe, conversion by Matt Weeks
 DATE CREATED: 27/5/94

 Harmonic sum pitch refinement function.

 pmin   pitch search range minimum
 pmax	pitch search range maximum
 step   pitch search step size
 model	current pitch estimate in model.Wo

 model 	refined pitch estimate in model.Wo

\*---------------------------------------------------------------------------*/
fn hs_pitch_refinement(model: &mut MODEL, Sw: &[COMP], pmin: f32, pmax: f32, pstep: f32) {
    //  Initialisation

    model.L = (PI / model.Wo as f64) as usize; //  use initial pitch est. for L
    let mut Wom = model.Wo; // Wo that maximises E
    let mut Em = 0.0; // mamimum energy
    let r = TWO_PI / FFT_ENC as f32; // number of rads/bin
    let one_on_r = 1.0 / r;

    //  Determine harmonic sum for a range of Wo values
    let mut p = pmin; // current pitch
    while p <= pmax {
        let mut E = 0.0; //energy for current pitch
        let Wo = TWO_PI / p; // current "test" fundamental freq.

        //  Sum harmonic magnitudes
        for m in 1..model.L + 1 {
            // bin for current harmonic centre
            let b = (m as f32 * Wo * one_on_r + 0.5) as usize;
            E += Sw[b].r * Sw[b].r + Sw[b].i * Sw[b].i;
        }
        //  Compare to see if this is a maximum

        if E > Em {
            Em = E;
            Wom = Wo;
        }
        p += pstep;
    }

    model.Wo = Wom;
}
