use crate::*;
const fdmdv_os_filter: [f32; 48] = [
    -0.0008215855034550382,
    -0.0007833023901802921,
    0.001075563790768233,
    0.001199092367787555,
    -0.001765309502928316,
    -0.002055372115328064,
    0.002986877604154257,
    0.003462567920638414,
    -0.004856570111126334,
    -0.005563143845031497,
    0.007533613299748122,
    0.008563932468880897,
    -0.01126857129039911,
    -0.01280782411693687,
    0.01651443896361847,
    0.01894875110322284,
    -0.02421604439474981,
    -0.02845107338464062,
    0.03672973563400258,
    0.04542046150312214,
    -0.06189165826716491,
    -0.08721876380763803,
    0.1496157094199961,
    0.4497962274137046,
    0.4497962274137046,
    0.1496157094199961,
    -0.08721876380763803,
    -0.0618916582671649,
    0.04542046150312216,
    0.03672973563400257,
    -0.02845107338464062,
    -0.02421604439474984,
    0.01894875110322284,
    0.01651443896361848,
    -0.01280782411693687,
    -0.0112685712903991,
    0.008563932468880899,
    0.007533613299748123,
    -0.005563143845031501,
    -0.004856570111126346,
    0.003462567920638419,
    0.002986877604154259,
    -0.002055372115328063,
    -0.001765309502928318,
    0.001199092367787557,
    0.001075563790768233,
    -0.0007833023901802925,
    -0.0008215855034550383,
];

const nlp_fir: [f32; 48] = [
    -1.0818124e-03,
    -1.1008344e-03,
    -9.2768838e-04,
    -4.2289438e-04,
    5.5034190e-04,
    2.0029849e-03,
    3.7058509e-03,
    5.1449415e-03,
    5.5924666e-03,
    4.3036754e-03,
    8.0284511e-04,
    -4.8204610e-03,
    -1.1705810e-02,
    -1.8199275e-02,
    -2.2065282e-02,
    -2.0920610e-02,
    -1.2808831e-02,
    3.2204775e-03,
    2.6683811e-02,
    5.5520624e-02,
    8.6305944e-02,
    1.1480192e-01,
    1.3674206e-01,
    1.4867556e-01,
    1.4867556e-01,
    1.3674206e-01,
    1.1480192e-01,
    8.6305944e-02,
    5.5520624e-02,
    2.6683811e-02,
    3.2204775e-03,
    -1.2808831e-02,
    -2.0920610e-02,
    -2.2065282e-02,
    -1.8199275e-02,
    -1.1705810e-02,
    -4.8204610e-03,
    8.0284511e-04,
    4.3036754e-03,
    5.5924666e-03,
    5.1449415e-03,
    3.7058509e-03,
    2.0029849e-03,
    5.5034190e-04,
    -4.2289438e-04,
    -9.2768838e-04,
    -1.1008344e-03,
    -1.0818124e-03,
];

//const PMAX_M: i32 = 320; /* maximum NLP analysis window size     */
const COEFF: f32 = 0.95; /* notch filter parameter               */
const PE_FFT_SIZE: i32 = 512; /* DFT size for pitch estimation        */
const DEC: i32 = 5; /* decimation factor                    */
const SAMPLE_RATE: i32 = 8000;
//const PI: f32 = 3.141592654; /* mathematical constant                */
//const T: f32 = 0.1; /* threshold for local minima candidate */
//const F0_MAX: i32 = 500;
const CNLP: f32 = 0.3; /* post processor constant              */
const NLP_NTAP: i32 = 48; /* Decimation LPF order */
const FDMDV_OS: i32 = 2; /* oversampling rate                   */
pub const FDMDV_OS_TAPS_16K: i32 = 48; /* number of OS filter taps at 16kHz   */
//const FDMDV_OS_TAPS_8K: i32 = FDMDV_OS_TAPS_16K / FDMDV_OS; /* number of OS filter taps at 8kHz    */
/*---------------------------------------------------------------------------*\

  post_process_sub_multiples()

  Given the global maximma of Fw[] we search integer submultiples for
  local maxima.  If local maxima exist and they are above an
  experimentally derived threshold (OK a magic number I pulled out of
  the air) we choose the submultiple as the F0 estimate.

  The rational for this is that the lowest frequency peak of Fw[]
  should be F0, as Fw[] can be considered the autocorrelation function
  of Sw[] (the speech spectrum).  However sometimes due to phase
  effects the lowest frequency maxima may not be the global maxima.

  This works OK in practice and favours low F0 values in the presence
  of background noise which means the sinusoidal codec does an OK job
  of synthesising the background noise.  High F0 in background noise
  tends to sound more periodic introducing annoying artifacts.

\*---------------------------------------------------------------------------*/

fn post_process_sub_multiples(
    Fw: &[COMP],
    _pmin: i32,
    pmax: i32,
    gmax: f32,
    gmax_bin: i32,
    prev_f0: &f32,
) -> f32 {
    /* post process estimate by searching submultiples */

    let mut mult = 2;
    let min_bin = PE_FFT_SIZE * DEC / pmax;
    let mut cmax_bin = gmax_bin;
    let prev_f0_bin = (*prev_f0) * ((PE_FFT_SIZE * DEC) as f32) / (SAMPLE_RATE as f32);

    while gmax_bin / mult >= min_bin {
        let b = gmax_bin / mult; /* determine search interval */
        let mut bmin = (0.8 * (b as f32)) as usize;
        let bmax = (1.2 * (b as f32)) as usize;
        if bmin < min_bin as usize {
            bmin = min_bin as usize;
        }
        /* lower threshold to favour previous frames pitch estimate,
        this is a form of pitch tracking */

        let thresh = if (prev_f0_bin > bmin as f32) && (prev_f0_bin < bmax as f32) {
            CNLP * 0.5 * gmax
        } else {
            CNLP * gmax
        };
        let mut lmax = 0.0;
        let mut lmax_bin = bmin;
        for b in bmin..bmax {
            /* look for maximum in interval */
            if Fw[b].r > lmax {
                lmax = Fw[b].r;
                lmax_bin = b;
            }
        }

        if lmax > thresh {
            if (lmax > Fw[lmax_bin - 1].r) && (lmax > Fw[lmax_bin + 1].r) {
                cmax_bin = lmax_bin as i32;
            }
        }
        mult += 1;
    }

    cmax_bin as f32 * SAMPLE_RATE as f32 / (PE_FFT_SIZE as f32 * DEC as f32)
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: fdmdv_16_to_8()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 9 May 2012

  Changes the sample rate of a signal from 16 to 8 kHz.

  n is the number of samples at the 8 kHz rate, there are FDMDV_OS*n
  samples at the 48 kHz rate.  As above however a memory of
  FDMDV_OS_TAPS samples is reqd for in16k[] (see t16_8.c unit test as example).

  Low pass filter the 16 kHz signal at 4 kHz using the same filter as
  the upsampler, then just output every FDMDV_OS-th filtered sample.

  Note: this function copied from fdmdv.c, included in nlp.c as a convenience
  to avoid linking with another source file.

\*---------------------------------------------------------------------------*/

fn fdmdv_16_to_8(out8k: &mut [f32], in16k: &mut [f32], in16koff: usize, n: i32) {
    let mut i = 0;
    for k in 0..n as usize {
        let mut acc = 0.0;
        for j in 0..FDMDV_OS_TAPS_16K as usize {
            acc += fdmdv_os_filter[j] * in16k[in16koff + i - j];
        }
        out8k[k] = acc;
        i += FDMDV_OS as usize;
    }

    /* update filter memory */

    for i in -FDMDV_OS_TAPS_16K..0 {
        in16k[(in16koff as i32 + i) as usize] =
            in16k[(in16koff as isize + i as isize + n as isize * FDMDV_OS as isize) as usize];
    }
}
/*---------------------------------------------------------------------------*\

  nlp()

  Determines the pitch in samples using the Non Linear Pitch (NLP)
  algorithm [1]. Returns the fundamental in Hz.  Note that the actual
  pitch estimate is for the centre of the M sample Sn[] vector, not
  the current N sample input vector.  This is (I think) a delay of 2.5
  frames with N=80 samples.  You should align further analysis using
  this pitch estimate to be centred on the middle of Sn[].

  Two post processors have been tried, the MBE version (as discussed
  in [1]), and a post processor that checks sub-multiples.  Both
  suffer occasional gross pitch errors (i.e. neither are perfect).  In
  the presence of background noise the sub-multiple algorithm tends
  towards low F0 which leads to better sounding background noise than
  the MBE post processor.

  A good way to test and develop the NLP pitch estimator is using the
  tnlp (codec2/unittest) and the codec2/octave/plnlp.m Octave script.

  A pitch tracker searching a few frames forward and backward in time
  would be a useful addition.

  References:

    [1] http://rowetel.com/downloads/1997_rowe_phd_thesis.pdf Chapter 4

\*---------------------------------------------------------------------------*/

pub fn nlp(
    nlp: &mut NLP,
    Sn: &[f32],        /* input speech vector                                */
    mut n: usize,      /* frames shift (no. new samples in Sn[])             */
    pitch: &mut f32,   /* estimated pitch period in samples at current Fs    */
    _Sw: &[COMP],      /* Freq domain version of Sn[]                        */
    _W: &[f32],        /* Freq domain window                                 */
    prev_f0: &mut f32, /* previous pitch f0 in Hz, memory for pitch tracking */
) -> f32 {
    //    notch: f32;		    /* current notch filter output          */
    //    let mut Fw[PE_FFT_SIZE];	    /* DFT of squared signal (input/output) */
    //    PROFILE_VAR(start, tnotch, filter, peakpick, window, fft, magsq, shiftmem);
    //    assert(nlp_state != NULL);
    //    nlp = (NLP*)nlp_state;
    let mut m = nlp.m;

    /* Square, notch filter at DC, and LP filter vector */

    /* If running at 16 kHz decimate to 8 kHz, as NLP ws designed for
    Fs = 8kHz. The decimating filter introduces about 3ms of delay,
    that shouldn't be a problem as pitch changes slowly. */

    if nlp.Fs == 8000 {
        /* Square latest input samples */

        for i in m - n..m {
            nlp.sq[i] = Sn[i] * Sn[i];
        }
    } else {
        //assert(nlp.Fs == 16000);

        /* re-sample at 8 KHz */

        for i in 0..n {
            nlp.Sn16k[FDMDV_OS_TAPS_16K as usize + i] = Sn[m + i - n];
        }

        m /= 2;
        n /= 2;

        let mut Sn8k = vec![0.0; n];
        fdmdv_16_to_8(
            &mut Sn8k,
            &mut nlp.Sn16k,
            FDMDV_OS_TAPS_16K as usize,
            n as i32,
        );

        /* Square latest input samples */
        let mut j = 0;
        for i in m - n..m {
            nlp.sq[i] = Sn8k[j] * Sn8k[j];
            j += 1;
        }
        //        assert(j <= n);
    }
    //fprintf(stderr, "n: %d m: %d\n", n, m);

    //    PROFILE_SAMPLE(start);
    for i in m - n..m {
        /* notch filter at DC */
        let mut notch = nlp.sq[i] - nlp.mem_x;
        notch += COEFF * nlp.mem_y;
        nlp.mem_x = nlp.sq[i];
        nlp.mem_y = notch;
        nlp.sq[i] = notch + 1.0; /* With 0 input vectors to codec,
                                 kiss_fft() would take a long
                                 time to execute when running in
                                 real time.  Problem was traced
                                 to kiss_fft function call in
                                 this function. Adding this small
                                 constant fixed problem.  Not
                                 exactly sure why. */
    }

    //    PROFILE_SAMPLE_AND_LOG(tnotch, start, "      square and notch");

    for i in m - n..m {
        /* FIR filter vector */

        for j in 0..NLP_NTAP as usize - 1 {
            nlp.mem_fir[j] = nlp.mem_fir[j + 1];
        }
        nlp.mem_fir[NLP_NTAP as usize - 1] = nlp.sq[i];

        nlp.sq[i] = 0.0;
        for j in 0..NLP_NTAP as usize {
            nlp.sq[i] += nlp.mem_fir[j] * nlp_fir[j];
        }
    }

    //    PROFILE_SAMPLE_AND_LOG(filter, tnotch, "      filter");

    /* Decimate and DFT */

    let mut Fw = [COMP::new(); PE_FFT_SIZE as usize];
    for i in 0..m / DEC as usize {
        Fw[i].r = nlp.sq[i * DEC as usize] * nlp.w[i];
    }
    //    PROFILE_SAMPLE_AND_LOG(window, filter, "      window");
    //    #ifdef DUMP
    //    dump_dec(Fw);
    //    #endif

    // FIXME: check if this can be converted to a real fft
    // since all imag inputs are 0
    crate::codec2_fft_inplace(&nlp.fft_cfg, &mut Fw);
    //    PROFILE_SAMPLE_AND_LOG(fft, window, "      fft");

    for i in 0..PE_FFT_SIZE as usize {
        Fw[i].r = Fw[i].r * Fw[i].r + Fw[i].i * Fw[i].i;
    }

    //    PROFILE_SAMPLE_AND_LOG(magsq, fft, "      mag sq");
    //    #ifdef DUMP
    //    dump_sq(m, nlp.sq);
    //    dump_Fw(Fw);
    //    #endif

    /* todo: express everything in f0, as pitch in samples is dep on Fs */

    let pmin = (SAMPLE_RATE as f32 * P_MIN_S).floor() as i32;
    let pmax = (SAMPLE_RATE as f32 * P_MAX_S).floor() as i32;

    /* find global peak */

    let mut gmax = 0.0;
    let mut gmax_bin = PE_FFT_SIZE * DEC / pmax;
    for i in (PE_FFT_SIZE * DEC / pmax) as usize..(PE_FFT_SIZE * DEC / pmin) as usize + 1 {
        if Fw[i].r > gmax {
            gmax = Fw[i].r;
            gmax_bin = i as i32;
        }
    }

    //    PROFILE_SAMPLE_AND_LOG(peakpick, magsq, "      peak pick");

    let best_f0 = post_process_sub_multiples(&Fw, pmin, pmax, gmax, gmax_bin, prev_f0);

    //    PROFILE_SAMPLE_AND_LOG(shiftmem, peakpick,  "      post process");

    /* Shift samples in buffer to make room for new samples */

    for i in 0..m - n {
        nlp.sq[i] = nlp.sq[i + n];
    }

    /* return pitch period in samples and F0 estimate */

    *pitch = nlp.Fs as f32 / best_f0;

    //    PROFILE_SAMPLE_AND_LOG2(shiftmem,  "      shift mem");

    //    PROFILE_SAMPLE_AND_LOG2(start,  "      nlp int");

    *prev_f0 = best_f0;

    best_f0
}
