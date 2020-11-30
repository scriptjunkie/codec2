/*
Copyright (c) 2003-2010, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

use crate::*;

/*
 Explanation of functions dealing with complex math:

  C_MUL(m,a,b)         : m = a*b
  C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
  C_SUB( res, a,b)     : res = a - b
  C_SUBFROM( res , a)  : res -= a
  C_ADDTO( res , a)    : res += a
* */

fn C_MUL(a: &kiss_fft_cpx, b: &kiss_fft_cpx) -> kiss_fft_cpx {
    kiss_fft_cpx {
        r: a.r * b.r - a.i * b.i,
        i: a.r * b.i + a.i * b.r,
    }
}
//#   define C_FIXDIV(c,div) /* NOOP */
fn C_MULBYSCALAR(c: &mut kiss_fft_cpx, s: f32) {
    c.r *= s;
    c.i *= s;
}
/*#ifndef CHECK_OVERFLOW_OP
#  define CHECK_OVERFLOW_OP(a,op,b) /* noop */
#endif*/

fn C_ADD(a: &kiss_fft_cpx, b: &kiss_fft_cpx) -> kiss_fft_cpx {
    kiss_fft_cpx {
        r: a.r + b.r,
        i: a.i + b.i,
    }
}
fn C_SUB(a: &kiss_fft_cpx, b: &kiss_fft_cpx) -> kiss_fft_cpx {
    kiss_fft_cpx {
        r: a.r - b.r,
        i: a.i - b.i,
    }
}
fn C_ADDTO(res: &mut kiss_fft_cpx, a: &kiss_fft_cpx) {
    res.r += a.r;
    res.i += a.i;
}

fn HALF_OF(x: f32) -> f32 {
    x / 2.0
}

fn kf_bfly2(Fout: &mut [kiss_fft_cpx], fstride: usize, st: &kiss_fft_cfg, mut m: usize) {
    let mut tw1 = &st.twiddles[..];
    let (mut Fout, mut Fout2) = Fout.split_at_mut(m);
    while m > 0 {
        //C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        let t = C_MUL(&Fout2[0], &tw1[0]);
        tw1 = &tw1[fstride..];
        Fout2[0] = C_SUB(&Fout[0], &t);
        C_ADDTO(&mut Fout[0], &t);
        Fout2 = &mut Fout2[1..];
        Fout = &mut Fout[1..];
        m = m - 1;
    }
}

fn kf_bfly4(mut Fout: &mut [kiss_fft_cpx], fstride: usize, st: &kiss_fft_cfg, m: usize) {
    let mut scratch = [kiss_fft_cpx::new(); 6];
    let mut k = m;
    let m2 = 2 * m;
    let m3 = 3 * m;

    let mut tw1 = &st.twiddles[..];
    let mut tw2 = &st.twiddles[..];
    let mut tw3 = &st.twiddles[..];

    loop {
        //        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        scratch[0] = C_MUL(&Fout[m], &tw1[0]);
        scratch[1] = C_MUL(&Fout[m2], &tw2[0]);
        scratch[2] = C_MUL(&Fout[m3], &tw3[0]);

        scratch[5] = C_SUB(&Fout[0], &scratch[1]);
        C_ADDTO(&mut Fout[0], &scratch[1]);
        scratch[3] = C_ADD(&scratch[0], &scratch[2]);
        scratch[4] = C_SUB(&scratch[0], &scratch[2]);
        Fout[m2] = C_SUB(&Fout[0], &scratch[3]);
        tw1 = &tw1[fstride..];
        tw2 = &tw2[fstride * 2..];
        tw3 = &tw3[fstride * 3..];
        C_ADDTO(&mut Fout[0], &scratch[3]);

        if st.inverse != 0 {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        } else {
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        Fout = &mut Fout[1..];
        k -= 1;
        if k == 0 {
            break;
        }
    }
}

fn kf_bfly3(mut Fout: &mut [kiss_fft_cpx], fstride: usize, st: &kiss_fft_cfg, m: usize) {
    let mut k = m;
    let m2 = 2 * m;
    let mut scratch = [kiss_fft_cpx::new(); 5];
    let epi3 = st.twiddles[fstride * m];

    let mut tw1 = 0;
    let mut tw2 = 0;

    while k > 0 {
        //         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

        scratch[1] = C_MUL(&Fout[m], &st.twiddles[tw1]);
        scratch[2] = C_MUL(&Fout[m2], &st.twiddles[tw2]);
        scratch[3] = C_ADD(&scratch[1], &scratch[2]);
        scratch[0] = C_SUB(&scratch[1], &scratch[2]);
        tw1 += fstride;
        tw2 += fstride * 2;

        Fout[m].r = Fout[0].r - (scratch[3].r / 2.0);
        Fout[m].i = Fout[0].i - (scratch[3].i / 2.0);

        C_MULBYSCALAR(&mut scratch[0], epi3.i);

        C_ADDTO(&mut Fout[0], &scratch[3]);

        Fout[m2].r = Fout[m].r + scratch[0].i;
        Fout[m2].i = Fout[m].i - scratch[0].r;

        Fout[m].r -= scratch[0].i;
        Fout[m].i += scratch[0].r;

        Fout = &mut Fout[1..];
        k -= 1;
    }
}

fn kf_bfly5(Fout: &mut [kiss_fft_cpx], fstride: usize, st: &kiss_fft_cfg, m: usize) {
    let mut scratch = [kiss_fft_cpx::new(); 13];
    let ya = st.twiddles[fstride * m];
    let yb = st.twiddles[fstride * 2 * m];

    let (mut Fout0, Fout1) = Fout.split_at_mut(m);
    let (mut Fout1, Fout2) = Fout1.split_at_mut(m);
    let (mut Fout2, Fout3) = Fout2.split_at_mut(m);
    let (mut Fout3, mut Fout4) = Fout3.split_at_mut(m);

    let tw = &st.twiddles[..];
    for u in 0..m {
        //        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = Fout0[0];

        scratch[1] = C_MUL(&Fout1[0], &tw[u * fstride]);
        scratch[2] = C_MUL(&Fout2[0], &tw[2 * u * fstride]);
        scratch[3] = C_MUL(&Fout3[0], &tw[3 * u * fstride]);
        scratch[4] = C_MUL(&Fout4[0], &tw[4 * u * fstride]);

        scratch[7] = C_ADD(&scratch[1], &scratch[4]);
        scratch[10] = C_SUB(&scratch[1], &scratch[4]);
        scratch[8] = C_ADD(&scratch[2], &scratch[3]);
        scratch[9] = C_SUB(&scratch[2], &scratch[3]);

        Fout0[0].r += scratch[7].r + scratch[8].r;
        Fout0[0].i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + (scratch[7].r * ya.r) + (scratch[8].r * yb.r);
        scratch[5].i = scratch[0].i + (scratch[7].i * ya.r) + (scratch[8].i * yb.r);

        scratch[6].r = (scratch[10].i * ya.i) + (scratch[9].i * yb.i);
        scratch[6].i = -(scratch[10].r * ya.i) - (scratch[9].r * yb.i);

        Fout1[0] = C_SUB(&scratch[5], &scratch[6]);
        Fout4[0] = C_ADD(&scratch[5], &scratch[6]);

        scratch[11].r = scratch[0].r + (scratch[7].r * yb.r) + (scratch[8].r * ya.r);
        scratch[11].i = scratch[0].i + (scratch[7].i * yb.r) + (scratch[8].i * ya.r);
        scratch[12].r = -(scratch[10].i * yb.i) + (scratch[9].i * ya.i);
        scratch[12].i = (scratch[10].r * yb.i) - (scratch[9].r * ya.i);

        Fout2[0] = C_ADD(&scratch[11], &scratch[12]);
        Fout3[0] = C_SUB(&scratch[11], &scratch[12]);

        Fout0 = &mut Fout0[1..];
        Fout1 = &mut Fout1[1..];
        Fout2 = &mut Fout2[1..];
        Fout3 = &mut Fout3[1..];
        Fout4 = &mut Fout4[1..];
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
fn kf_bfly_generic(
    Fout: &mut [kiss_fft_cpx],
    fstride: usize,
    st: &kiss_fft_cfg,
    m: usize,
    p: usize,
) {
    let Norig = st.nfft;

    let mut scratch = vec![kiss_fft_cpx::new(); p];

    for u in 0..m {
        let mut k = u;
        for q1 in 0..p {
            scratch[q1] = Fout[k];
            //            C_FIXDIV(scratch[q1],p);
            k += m;
        }

        k = u;
        for _q1 in 0..p {
            let mut twidx = 0;
            Fout[k] = scratch[0];
            for q in 1..p {
                twidx += fstride * k;
                if twidx >= Norig as usize {
                    twidx -= Norig as usize;
                }
                let t = C_MUL(&scratch[q], &st.twiddles[twidx]);
                C_ADDTO(&mut Fout[k], &t);
            }
            k += m;
        }
    }
}

fn kf_work(
    Fout: &mut [kiss_fft_cpx],
    mut f: &[kiss_fft_cpx],
    fstride: usize,
    in_stride: i32,
    factors: &[i32],
    st: &kiss_fft_cfg,
) {
    let p = factors[0] as usize; /* the radix  */
    let m = factors[1] as usize; /* stage's fft length/p */
    if m == 1 {
        for fout_idx in 0..p * m {
            Fout[fout_idx] = f[0];
            f = &f[(fstride * in_stride as usize).min(f.len())..];
        }
    } else {
        for fout_idx in (0..p * m).step_by(m) {
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work(
                &mut Fout[fout_idx..],
                f,
                fstride * p,
                in_stride,
                &factors[2..],
                st,
            );
            f = &f[(fstride * in_stride as usize).min(f.len())..];
        }
    }

    // recombine the p smaller DFTs
    match p {
        2 => {
            kf_bfly2(Fout, fstride, st, m);
        }
        3 => {
            kf_bfly3(Fout, fstride, st, m);
        }
        4 => {
            kf_bfly4(Fout, fstride, st, m);
        }
        5 => {
            kf_bfly5(Fout, fstride, st, m);
        }
        _ => {
            kf_bfly_generic(Fout, fstride, st, m, p);
        }
    }
}

fn kiss_fft_stride(
    st: &kiss_fft_cfg,
    fin: &[kiss_fft_cpx],
    fout: &mut [kiss_fft_cpx],
    in_stride: i32,
) {
    kf_work(fout, fin, 1, in_stride, &st.factors[..], st);
}

pub fn kiss_fft(cfg: &kiss_fft_cfg, fin: &[kiss_fft_cpx], fout: &mut [kiss_fft_cpx]) {
    kiss_fft_stride(cfg, fin, fout, 1);
}

pub fn kiss_fftr(st: &mut kiss_fftr_cfg, timedata: &[f32], freqdata: &mut [kiss_fft_cpx]) {
    /* input buffer timedata is stored row-wise */
    //    assert(st.substate.inverse==0);

    let ncfft = st.substate.nfft as usize;

    /*perform the parallel fft of two real signals packed in real,imag*/
    //could transmute this but UGH so we copy
    let mut timedatac = vec![kiss_fft_cpx::new(); timedata.len() / 2];
    for i in 0..timedatac.len() / 2 {
        timedatac[i] = kiss_fft_cpx {
            r: timedata[i * 2],
            i: timedata[i * 2 + 1],
        };
    }
    kiss_fft(&st.substate, &timedatac, &mut st.tmpbuf);

    /* The real part of the DC element of the frequency spectrum in st.tmpbuf
     * contains the sum of the even-numbered elements of the input time sequence
     * The imag part is the sum of the odd-numbered elements
     *
     * The sum of tdc.r and tdc.i is the sum of the input time sequence.
     *      yielding DC of input time sequence
     * The difference of tdc.r - tdc.i is the sum of the input (dot product) [1,-1,1,-1...
     *      yielding Nyquist bin of input time sequence
     */
    let mut tdc = kiss_fft_cpx::new();
    tdc.r = st.tmpbuf[0].r;
    tdc.i = st.tmpbuf[0].i;
    //    C_FIXDIV(tdc,2);
    //    CHECK_OVERFLOW_OP(tdc.r ,+, tdc.i);
    //    CHECK_OVERFLOW_OP(tdc.r ,-, tdc.i);
    freqdata[0].r = tdc.r + tdc.i;
    freqdata[ncfft].r = tdc.r - tdc.i;
    freqdata[ncfft].i = 0.0;
    freqdata[0].i = 0.0;

    for k in 1..ncfft / 2 + 1 {
        let fpk = st.tmpbuf[k];
        let mut fpnk = kiss_fft_cpx::new();
        fpnk.r = st.tmpbuf[ncfft - k].r;
        fpnk.i = 0.0 - st.tmpbuf[ncfft - k].i;
        //        C_FIXDIV(fpk,2);
        //        C_FIXDIV(fpnk,2);

        let f1k = C_ADD(&fpk, &fpnk);
        let f2k = C_SUB(&fpk, &fpnk);
        let tw = C_MUL(&f2k, &st.super_twiddles[k - 1]);

        freqdata[k].r = HALF_OF(f1k.r + tw.r);
        freqdata[k].i = HALF_OF(f1k.i + tw.i);
        freqdata[ncfft - k].r = HALF_OF(f1k.r - tw.r);
        freqdata[ncfft - k].i = HALF_OF(tw.i - f1k.i);
    }
}

pub fn kiss_fftri(st: &mut kiss_fftr_cfg, freqdata: &[kiss_fft_cpx], timedata: &mut [f32]) {
    /* input buffer timedata is stored row-wise */
    //   assert(st.substate.inverse == 1);

    let ncfft = st.substate.nfft as usize;

    st.tmpbuf[0].r = freqdata[0].r + freqdata[ncfft].r;
    st.tmpbuf[0].i = freqdata[0].r - freqdata[ncfft].r;
    //    C_FIXDIV(st.tmpbuf[0],2);

    for k in 1..ncfft / 2 + 1 {
        let fk = freqdata[k];
        let mut fnkc = kiss_fft_cpx::new();
        fnkc.r = freqdata[ncfft - k].r;
        fnkc.i = 0.0 - freqdata[ncfft - k].i;
        //        C_FIXDIV( fk , 2 );
        //        C_FIXDIV( fnkc , 2 );

        let fek = C_ADD(&fk, &fnkc);
        let tmp = C_SUB(&fk, &fnkc);
        let fok = C_MUL(&tmp, &st.super_twiddles[k - 1]);
        st.tmpbuf[k] = C_ADD(&fek, &fok);
        st.tmpbuf[ncfft - k] = C_SUB(&fek, &fok);
        st.tmpbuf[ncfft - k].i *= -1.0;
    }
    //could transmute this but UGH so we copy
    let mut timedatac = vec![kiss_fft_cpx::new(); timedata.len() / 2];
    kiss_fft(&st.substate, &st.tmpbuf, &mut timedatac);
    for i in 0..timedatac.len() {
        timedata[i * 2] = timedatac[i].r;
        timedata[i * 2 + 1] = timedatac[i].i;
    }
}
