use crate::*;
const LSP_DELTA1: f32 = 0.01;
const E_MIN_DB: f32 = -10.0;
const E_MAX_DB: f32 = 40.0;
pub const E_BITS: i32 = 5;

pub fn lspd_bits(i: usize) -> i32 {
    lsp_cbd[i].log2m
}

/*---------------------------------------------------------------------------*\

  autocorrelate()

  Finds the first P autocorrelation values of an array of windowed speech
  samples Sn[].

\*---------------------------------------------------------------------------*/
pub fn autocorrelate(
    Sn: &[f32],     //  frame of Nsam windowed speech samples
    Rn: &mut [f32], //  array of P+1 autocorrelation coefficients
    Nsam: usize,    //  number of windowed samples to use
    order: usize,   //  order of LPC analysis
) {
    for j in 0..order + 1 {
        Rn[j] = 0.0;
        for i in 0..Nsam - j {
            Rn[j] += Sn[i] * Sn[i + j];
        }
    }
}

/*---------------------------------------------------------------------------*\

  levinson_durbin()

  Given P+1 autocorrelation coefficients, finds P Linear Prediction Coeff.
  (LPCs) where P is the order of the LPC all-pole model. The Levinson-Durbin
  algorithm is used, and is described in:

    J. Makhoul
    "Linear prediction, a tutorial review"
    Proceedings of the IEEE
    Vol-63, No. 4, April 1975

\*---------------------------------------------------------------------------*/
pub fn levinson_durbin(
    R: &[f32],        //  order+1 autocorrelation coeff
    lpcs: &mut [f32], //  order+1 LPC's
    order: usize,     //  order of the LPC analysis
) {
    let mut a = vec![vec![0.0; order + 1]; order + 1];
    let mut e = R[0]; //  Equation 38a, Makhoul

    for i in 1..order + 1 {
        let mut sum = 0.0;
        for j in 1..i {
            sum += a[i - 1][j] * R[i - j];
        }
        let mut k = -1.0 * (R[i] + sum) / e; //  Equation 38b, Makhoul
        if k.abs() > 1.0 {
            k = 0.0;
        }
        a[i][i] = k;

        for j in 1..i {
            a[i][j] = a[i - 1][j] + k * a[i - 1][i - j]; //  Equation 38c, Makhoul
        }
        e *= 1.0 - k * k; //  Equation 38d, Makhoul
    }

    for i in 1..order + 1 {
        lpcs[i] = a[order][i];
    }
    lpcs[0] = 1.0;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: cheb_poly_eva()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 24/2/93

  This function evalutes a series of chebyshev polynomials

  FIXME: performing memory allocation at run time is very inefficient,
  replace with stack variables of MAX_P size.

\*---------------------------------------------------------------------------*/
fn cheb_poly_eva(
    coef: &[f32], //   coefficients of the polynomial to be evaluated
    x: f32,       //  the point where polynomial is to be evaluated
    order: usize, //  order of the polynomial
) -> f32 {
    let mut T = vec![0.0; order / 2 + 1];

    //  Initialise offsets
    T[0] = 1.0;
    T[1] = x;
    let mut toff = 0;
    let mut uoff = 1;
    let mut voff = 2;
    //  Evaluate chebyshev series formulation using iterative approach

    for _i in 2..order / 2 + 1 {
        T[voff] = (2.0 * x) * T[uoff] - T[toff]; //  T[i] = 2*x*T[i-1] - T[i-2]
        toff += 1;
        uoff += 1;
        voff += 1;
    }
    let mut sum = 0.0; //  initialise sum to zero
    toff = 0; //  reset offset

    //  Evaluate polynomial and return value also free memory space

    for i in 0..order / 2 + 1 {
        sum += coef[(order / 2) - i] * T[toff];
        toff += 1;
    }
    return sum;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: lpc_to_lsp()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 24/2/93

  This function converts LPC coefficients to LSP coefficients.

\*---------------------------------------------------------------------------*/
pub fn lpc_to_lsp(
    a: &[f32],        //   float *a      lpc coefficients
    order: usize,     //   int order     order of LPC coefficients (10)
    freq: &mut [f32], //   float *freq   LSP frequencies in radians
    nb: i32,          //   int nb        number of sub-intervals (4)
    delta: f32,       //   float delta   grid spacing interval (0.02)
) -> i32 {
    let mut roots = 0; //  number of roots found
    let mut Q = vec![0.0; order + 1];
    let mut P = vec![0.0; order + 1];

    let m = order / 2; //  order of P'(z) & Q'(z) polynimials

    //  Allocate memory space for polynomials

    //  determine P'(z)'s and Q'(z)'s coefficients where
    // P'(z) = P(z)/(1 + z^(-1)) and Q'(z) = Q(z)/(1-z^(-1))
    let mut pxoff = 0; //offsets of respective P'(z) & Q'(z)
    let mut qxoff = 0;
    let mut poff = pxoff;
    let mut qoff = qxoff;
    P[pxoff] = 1.0;
    pxoff += 1;
    Q[qxoff] = 1.0;
    qxoff += 1;
    for i in 1..m + 1 {
        P[pxoff] = a[i] + a[order + 1 - i] - P[poff];
        pxoff += 1;
        poff += 1;
        Q[qxoff] = a[i] - a[order + 1 - i] + Q[qoff];
        qxoff += 1;
        qoff += 1;
    }
    pxoff = 0;
    qxoff = 0;
    for _i in 0..m {
        P[pxoff] = 2.0 * P[pxoff];
        Q[qxoff] = 2.0 * Q[qxoff];
        pxoff += 1;
        qxoff += 1;
    }

    //  Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
    //Keep alternating between the two polynomials as each zero is found

    let mut xr = 0.0; //  initialise xr to zero
    let mut xl = 1.0; //  start at point xl = 1

    for j in 0..order {
        //ptr used for cheb_poly_eval()
        let pt = if (j % 2) != 0 {
            //  determines whether P' or Q' is eval.
            &Q[..]
        } else {
            &P[..]
        };
        let mut psuml = cheb_poly_eva(pt, xl, order); //  evals poly. at xl
        let mut flag = true;
        while flag && (xr >= -1.0) {
            xr = xl - delta; //  interval spacing
            let psumr = cheb_poly_eva(pt, xr, order); //  poly(xl-delta_x)
            let temp_psumr = psumr;
            let temp_xr = xr;

            //  if no sign change increment xr and re-evaluate
            // poly(xr). Repeat til sign change.  if a sign change has
            // occurred the interval is bisected and then checked again
            // for a sign change which determines in which interval the
            // zero lies in.  If there is no sign change between poly(xm)
            // and poly(xl) set interval between xm and xr else set
            // interval between xl and xr and repeat till root is located
            // within the specified limits

            if ((psumr * psuml) < 0.0) || (psumr == 0.0) {
                roots += 1;

                let mut xm = 0.0;
                for _k in 0..nb + 1 {
                    xm = (xl + xr) / 2.0; //  bisect the interval
                    let psumm = cheb_poly_eva(pt, xm, order);
                    if psumm * psuml > 0. {
                        psuml = psumm;
                        xl = xm;
                    } else {
                        //psumr = psumm; //TODO: unused?
                        xr = xm;
                    }
                }

                //  once zero is found, reset initial interval to xr
                freq[j] = xm;
                xl = xm;
                flag = false; //  reset flag for next search
            } else {
                psuml = temp_psumr;
                xl = temp_xr;
            }
        }
    }

    //  convert from x domain to radians
    for i in 0..order {
        freq[i] = freq[i].acos();
    }

    return roots;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: encode_Wo()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Encodes Wo using a WO_LEVELS quantiser.

\*---------------------------------------------------------------------------*/
pub fn encode_Wo(c2const: &C2const, Wo: f32, bits: i32) -> i32 {
    let mut index = 1 << bits;
    let Wo_levels = index;
    let Wo_min = c2const.Wo_min;
    let Wo_max = c2const.Wo_max;
    let norm = (Wo - Wo_min) / (Wo_max - Wo_min);
    index = (Wo_levels as f32 * norm + 0.5).floor() as i32;
    if index < 0 {
        index = 0;
    }
    if index > (Wo_levels - 1) {
        index = Wo_levels - 1;
    }
    index
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: decode_Wo()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Decodes Wo using a WO_LEVELS quantiser.

\*---------------------------------------------------------------------------*/
pub fn decode_Wo(c2const: &C2const, index: i32, bits: i32) -> f32 {
    let Wo_min = c2const.Wo_min;
    let Wo_max = c2const.Wo_max;
    let Wo_levels = 1 << bits;

    let step = (Wo_max - Wo_min) / Wo_levels as f32;
    Wo_min + step * index as f32
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: speech_to_uq_lsps()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Analyse a windowed frame of time domain speech to determine LPCs
  which are the converted to LSPs for quantisation and transmission
  over the channel.

\*---------------------------------------------------------------------------*/
pub fn speech_to_uq_lsps(
    lsp: &mut [f32],
    ak: &mut [f32],
    Sn: &[f32],
    w: &[f32],
    m_pitch: usize,
    order: usize,
) -> f32 {
    let mut Wn = vec![0.0; m_pitch];
    let mut R = vec![0.0; order + 1];

    let mut e = 0.0;
    for i in 0..m_pitch {
        Wn[i] = Sn[i] * w[i];
        e += Wn[i] * Wn[i];
    }

    //  trap 0 energy case as LPC analysis will fail

    if e == 0.0 {
        for i in 0..order {
            lsp[i] = (PI as f32 / order as f32) * (i as f32);
        }
        return 0.0;
    }

    autocorrelate(&Wn, &mut R, m_pitch, order);
    levinson_durbin(&R, ak, order);

    let mut E = 0.0;
    for i in 0..order + 1 {
        E += ak[i] * R[i];
    }

    //  15 Hz BW expansion as I can't hear the difference and it may help
    //   help occasional fails in the LSP root finding.  Important to do this
    //   after energy calculation to avoid -ve energy values.

    for i in 0..order + 1 {
        ak[i] *= 0.994f32.powf(i as f32);
    }

    let roots = lpc_to_lsp(ak, order, lsp, 5, LSP_DELTA1);
    if roots != order as i32 {
        //  if root finding fails use some benign LSP values instead
        for i in 0..order {
            lsp[i] = (PI as f32 / order as f32) * (i as f32);
        }
    }

    return E;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: encode_energy()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Encodes LPC energy using an E_LEVELS quantiser.

\*---------------------------------------------------------------------------*/
pub fn encode_energy(e: f32, bits: i32) -> i32 {
    let e_levels = 1 << bits;

    let e = 10.0 * e.log10();
    let norm = (e - E_MIN_DB) / (E_MAX_DB - E_MIN_DB);
    let mut index = (e_levels as f32 * norm + 0.5).floor() as i32;
    if index < 0 {
        index = 0;
    }
    if index > (e_levels - 1) {
        index = e_levels - 1;
    }

    return index;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: decode_energy()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Decodes energy using a E_LEVELS quantiser.

\*---------------------------------------------------------------------------*/
pub fn decode_energy(index: i32, bits: i32) -> f32 {
    let e_min = E_MIN_DB;
    let e_max = E_MAX_DB;
    let e_levels = 1 << bits;

    let step = (e_max - e_min) / e_levels as f32;
    let e = e_min + step * (index as f32);
    10.0_f32.powf(e / 10.0)
}

/*---------------------------------------------------------------------------*\

  encode_lspds_scalar()

  Scalar/VQ LSP difference quantiser.

\*---------------------------------------------------------------------------*/
pub fn encode_lspds_scalar(indexes: &mut [i32], lsp: &[f32], order: usize) {
    let mut lsp_hz = vec![0.0; order];
    let mut lsp__hz = vec![0.0; order];
    let mut dlsp = vec![0.0; order];
    let mut dlsp_ = vec![0.0; order];
    let mut wt = vec![1.0; order];

    //  convert from radians to Hz so we can use human readable frequencies

    for i in 0..order {
        lsp_hz[i] = (4000.0 / PI as f32) * lsp[i];
    }
    wt[0] = 1.0;
    for i in 0..order {
        //  find difference from previous qunatised lsp

        if i != 0 {
            dlsp[i] = lsp_hz[i] - lsp__hz[i - 1];
        } else {
            dlsp[0] = lsp_hz[0];
        }
        let k = lsp_cbd[i].k as usize;
        let m = lsp_cbd[i].m as usize;
        let cb = lsp_cbd[i].cb;
        let mut se = 0.0;
        indexes[i] = quantise(cb, &mut dlsp[i..], &wt, k, m, &mut se) as i32;
        dlsp_[i] = cb[indexes[i] as usize * k];

        if i != 0 {
            lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
        } else {
            lsp__hz[0] = dlsp_[0];
        }
    }
}

/*---------------------------------------------------------------------------*\

  quantise

  Quantises vec by choosing the nearest vector in codebook cb, and
  returns the vector index.  The squared error of the quantised vector
  is added to se.

\*---------------------------------------------------------------------------*/
fn quantise(
    cb: &[f32],   //  cb[][K];  current VQ codebook
    vec: &[f32],  //  vec[];    vector to quantise
    w: &[f32],    //  w[];      weighting vector
    k: usize,     //  k;        dimension of vectors
    m: usize,     //  m;        size of codebook
    se: &mut f32, //  *se;      accumulated squared error
) -> i64 {
    let mut besti = 0; // best index so far
    let mut beste = 1E32; // best error so far
    for j in 0..m {
        let mut e = 0.0; // current error
        for i in 0..k {
            let diff = cb[j * k + i] - vec[i];
            e += diff * w[i] * diff * w[i];
        }
        if e < beste {
            beste = e;
            besti = j;
        }
    }

    *se += beste;
    besti as i64
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: interp_energy()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22 May 2012

  Interpolates centre 10ms sample of energy given two samples 20ms
  apart.

\*---------------------------------------------------------------------------*/
pub fn interp_energy(prev_e: f32, next_e: f32) -> f32 {
    //return powf(10.0, (log10f(prev_e) + log10f(next_e))/2.0);
    (prev_e * next_e).sqrt() //looks better is math. identical and faster math
}

/*---------------------------------------------------------------------------*\

   aks_to_M2()

   Transforms the linear prediction coefficients to spectral amplitude
   samples.  This function determines A(m) from the average energy per
   band using an FFT.

\*---------------------------------------------------------------------------*/
pub fn aks_to_M2(
    fftr_fwd_cfg: &mut codec2_fftr_cfg,
    ak: &[f32], //  LPC's
    order: usize,
    model: &mut MODEL, //  sinusoidal model parameters for this frame
    E: f32,            //  energy term
    snr: &mut f32,     //  signal to noise ratio for this frame in dB
    dump: i32,         //  true to dump sample to dump file
    sim_pf: i32,       //  true to simulate a post filter
    pf: i32,           //  true to enable actual LPC post filter
    bass_boost: i32,   //  enable LPC filter 0-1kHz 3dB boost
    beta: f32,
    gamma: f32,      //  LPC post filter parameters
    Aw: &mut [COMP], //  output power spectrum
) {
    //PROFILE_VAR(tstart, tfft, tpw, tpf);
    //
    //PROFILE_SAMPLE(tstart);

    let r = TWO_PI / FFT_ENC as f32; // no. rads/bin

    //  Determine DFT of A(exp(jw)) --------------------------------------------
    {
        let mut a = [0.0; FFT_ENC]; //  input to FFT for power spectrum

        for i in 0..order + 1 {
            a[i] = ak[i];
        }
        codec2_fftr(fftr_fwd_cfg, &a, Aw);
    }
    //  PROFILE_SAMPLE_AND_LOG(tfft, tstart, "      fft");

    //  Determine power spectrum P(w) = E/(A(exp(jw))^2 ------------------------

    let mut Pw = [0.0; FFT_ENC / 2];

    //#ifndef FDV_ARM_MATH
    for i in 0..FFT_ENC / 2 {
        Pw[i] = 1.0 / (Aw[i].r * Aw[i].r + Aw[i].i * Aw[i].i + 1E-6);
    }
    /* #else
      // this difference may seem strange, but the gcc for STM32F4 generates almost 5 times
      // faster code with the two loops: 1120 ms -> 242 ms
      // so please leave it as is or improve further
      // since this code is called 4 times it results in almost 4ms gain (21ms -> 17ms per audio frame decode @ 1300 )

      for i in 0..FFT_ENC/2
      {
          Pw[i] = Aw[i].r * Aw[i].r + Aw[i].i * Aw[i].i  + 1E-6;
      }
      for i in 0..FFT_ENC/2 {
          Pw[i] = 1.0/(Pw[i]);
      }
    #endif // */

    //  PROFILE_SAMPLE_AND_LOG(tpw, tfft, "      Pw");

    if pf != 0 {
        lpc_post_filter(
            fftr_fwd_cfg,
            &mut Pw,
            ak,
            order,
            dump,
            beta,
            gamma,
            bass_boost,
            E,
        );
    } else {
        for i in 0..FFT_ENC / 2 {
            Pw[i] *= E;
        }
    }

    //  PROFILE_SAMPLE_AND_LOG(tpf, tpw, "      LPC post filter");

    //  #ifdef DUMP
    //  if (dump)
    //      dump_Pw(Pw);
    //  #endif

    //  Determine magnitudes from P(w) ----------------------------------------

    //  when used just by decoder {A} might be all zeroes so init signal
    // and noise to prevent log(0) errors

    let mut signal = 1E-30;
    let mut noise = 1E-32;

    for m in 1..model.L as usize + 1 {
        //limits of current band
        let am = ((m as f32 - 0.5) * model.Wo / r + 0.5) as usize;
        let mut bm = ((m as f32 + 0.5) * model.Wo / r + 0.5) as usize;

        // FIXME: With arm_rfft_fast_f32 we have to use this
        // otherwise sometimes a to high bm is calculated
        // which causes trouble later in the calculation
        // chain
        // it seems for some reason model.Wo is calculated somewhat too high
        if bm > FFT_ENC / 2 {
            bm = FFT_ENC / 2;
        }
        let mut Em = 0.0_f32; // energy in band

        for i in am..bm {
            Em += Pw[i];
        }
        let mut Am = Em.sqrt(); //  spectral amplitude sample

        signal += model.A[m] * model.A[m];
        noise += (model.A[m] - Am) * (model.A[m] - Am);

        //  This code significantly improves perf of LPC model, in
        // particular when combined with phase0.  The LPC spectrum tends
        // to track just under the peaks of the spectral envelope, and
        // just above nulls.  This algorithm does the reverse to
        // compensate - raising the amplitudes of spectral peaks, while
        // attenuating the null.  This enhances the formants, and
        // supresses the energy between formants.

        if sim_pf != 0 {
            if Am > model.A[m] {
                Am *= 0.7;
            }
            if Am < model.A[m] {
                Am *= 1.4;
            }
        }
        model.A[m] = Am;
    }
    *snr = 10.0 * (signal / noise).log10();

    //  PROFILE_SAMPLE_AND_LOG2(tpf, "      rec");
}

fn cconj(a: COMP) -> COMP {
    COMP {
        r: a.r,
        i: 0.0 - a.i,
    }
}

/*---------------------------------------------------------------------------*\

  sample_phase()

  Samples phase at centre of each harmonic from and array of FFT_ENC
  DFT samples.

\*---------------------------------------------------------------------------*/
pub fn sample_phase(
    model: &MODEL,
    H: &mut [COMP],
    A: &[COMP], //  LPC analysis filter in freq domain
) {
    let r = TWO_PI / (FFT_ENC as f32);

    //  Sample phase at harmonics

    for m in 1..model.L as usize + 1 {
        let b = (m as f32 * model.Wo / r + 0.5) as usize;
        H[m] = cconj(A[b]); //  synth filter 1/A is opposite phase to analysis filter
    }
}

/*---------------------------------------------------------------------------*\

   phase_synth_zero_order()

   Synthesises phases based on SNR and a rule based approach.  No phase
   parameters are required apart from the SNR (which can be reduced to a
   1 bit V/UV decision per frame).

   The phase of each harmonic is modelled as the phase of a synthesis
   filter excited by an impulse.  In many Codec 2 modes the synthesis
   filter is a LPC filter. Unlike the first order model the position
   of the impulse is not transmitted, so we create an excitation pulse
   train using a rule based approach.

   Consider a pulse train with a pulse starting time n=0, with pulses
   repeated at a rate of Wo, the fundamental frequency.  A pulse train
   in the time domain is equivalent to harmonics in the frequency
   domain.  We can make an excitation pulse train using a sum of
   sinsusoids:

     for m in 1..L
       ex[n] = cos(m*Wo*n)

   Note: the Octave script ../octave/phase.m is an example of this if
   you would like to try making a pulse train.

   The phase of each excitation harmonic is:

     arg(E[m]) = mWo

   where E[m] are the complex excitation (freq domain) samples,
   arg(x), just returns the phase of a complex sample x.

   As we don't transmit the pulse position for this model, we need to
   synthesise it.  Now the excitation pulses occur at a rate of Wo.
   This means the phase of the first harmonic advances by N_SAMP samples
   over a synthesis frame of N_SAMP samples.  For example if Wo is pi/20
   (200 Hz), then over a 10ms frame (N_SAMP=80 samples), the phase of the
   first harmonic would advance (pi/20)*80 = 4*pi or two complete
   cycles.

   We generate the excitation phase of the fundamental (first
   harmonic):

     arg[E[1]] = Wo*N_SAMP;

   We then relate the phase of the m-th excitation harmonic to the
   phase of the fundamental as:

     arg(E[m]) = m*arg(E[1])

   This E[m] then gets passed through the LPC synthesis filter to
   determine the final harmonic phase.

   Comparing to speech synthesised using original phases:

   - Through headphones speech synthesised with this model is not as
     good. Through a loudspeaker it is very close to original phases.

   - If there are voicing errors, the speech can sound clicky or
     staticy.  If V speech is mistakenly declared UV, this model tends to
     synthesise impulses or clicks, as there is usually very little shift or
     dispersion through the LPC synthesis filter.

   - When combined with LPC amplitude modelling there is an additional
     drop in quality.  I am not sure why, theory is interformant energy
     is raised making any phase errors more obvious.

   NOTES:

     1/ This synthesis model is effectively the same as a simple LPC-10
     vocoders, and yet sounds much better.  Why? Conventional wisdom
     (AMBE, MELP) says mixed voicing is required for high quality
     speech.

     2/ I am pretty sure the Lincoln Lab sinusoidal coding guys (like xMBE
     also from MIT) first described this zero phase model, I need to look
     up the paper.

     3/ Note that this approach could cause some discontinuities in
     the phase at the edge of synthesis frames, as no attempt is made
     to make sure that the phase tracks are continuous (the excitation
     phases are continuous, but not the final phases after filtering
     by the LPC spectra).  Technically this is a bad thing.  However
     this may actually be a good thing, disturbing the phase tracks a
     bit.  More research needed, e.g. test a synthesis model that adds
     a small delta-W to make phase tracks line up for voiced
     harmonics.

\*---------------------------------------------------------------------------*/
pub fn phase_synth_zero_order(
    n_samp: i32,
    model: &mut MODEL,
    ex_phase: &mut f32, //  excitation phase of fundamental
    H: &[COMP],         //  L synthesis filter freq domain samples
) {
    let mut Ex = [COMP::new(); MAX_AMP + 1]; //  excitation samples
    let mut A_ = [COMP::new(); MAX_AMP + 1]; //  synthesised harmonic samples

    /*
       Update excitation fundamental phase track, this sets the position
       of each pitch pulse during voiced speech.  After much experiment
       I found that using just this frame's Wo improved quality for UV
       sounds compared to interpolating two frames Wo like this:

       *ex_phase += (*prev_Wo+model.Wo)*N_SAMP/2;
    */

    *ex_phase += (model.Wo) * n_samp as f32;
    *ex_phase -= TWO_PI * (*ex_phase / TWO_PI + 0.5).floor();

    for m in 1..model.L as usize + 1 {
        //  generate excitation

        if model.voiced != 0 {
            Ex[m].r = (*ex_phase * m as f32).cos();
            Ex[m].i = (*ex_phase * m as f32).sin();
        } else {
            //  When a few samples were tested I found that LPC filter
            //   phase is not needed in the unvoiced case, but no harm in
            //   keeping it.

            let phi = TWO_PI * codec2_rand() / CODEC2_RAND_MAX;
            Ex[m].r = phi.cos();
            Ex[m].i = phi.sin();
        }

        //  filter using LPC filter

        A_[m].r = H[m].r * Ex[m].r - H[m].i * Ex[m].i;
        A_[m].i = H[m].i * Ex[m].r + H[m].r * Ex[m].i;

        //  modify sinusoidal phase

        let new_phi = A_[m].i.atan2(A_[m].r + 1E-12);
        model.phi[m] = new_phi;
    }
}

const CODEC2_RAND_MAX: f32 = 32767.0;

//  todo: this should probably be in some states rather than a static
static next_rand: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(1);
fn codec2_rand() -> f32 {
    // next = next * 1103515245 + 12345 but allow overflow and allow thready accesses.
    let mut nextr = next_rand.load(std::sync::atomic::Ordering::Relaxed);
    nextr = nextr.overflowing_mul(1103515245).0.overflowing_add(12345).0;
    next_rand.store(nextr, std::sync::atomic::Ordering::Relaxed);
    ((nextr / 65536) % 32768) as f32
}

/*---------------------------------------------------------------------------*\

   lpc_post_filter()

   Applies a post filter to the LPC synthesis filter power spectrum
   Pw, which supresses the inter-formant energy.

   The algorithm is from p267 (Section 8.6) of "Digital Speech",
   edited by A.M. Kondoz, 1994 published by Wiley and Sons.  Chapter 8
   of this text is on the MBE vocoder, and this is a freq domain
   adaptation of post filtering commonly used in CELP.

   I used the Octave simulation lpcpf.m to get an understanding of the
   algorithm.

   Requires two more FFTs which is significantly more MIPs.  However
   it should be possible to implement this more efficiently in the
   time domain.  Just not sure how to handle relative time delays
   between the synthesis stage and updating these coeffs.  A smaller
   FFT size might also be accetable to save CPU.

   TODO:
   [ ] sync var names between Octave and C version
   [ ] doc gain normalisation
   [ ] I think the first FFT is not rqd as we do the same
       thing in aks_to_M2().

\*---------------------------------------------------------------------------*/
fn lpc_post_filter(
    fftr_fwd_cfg: &mut codec2_fftr_cfg,
    Pw: &mut [f32],
    ak: &[f32],
    order: usize,
    _dump: i32,
    beta: f32,
    gamma: f32,
    bass_boost: i32,
    E: f32,
) {
    let mut x = [0.0; FFT_ENC]; //  input to FFTs
    let mut Ww = [COMP::new(); FFT_ENC / 2 + 1]; //  weighting spectrum
    let mut Rw = [0.0; FFT_ENC / 2 + 1]; //  R = WA
                                         //    PROFILE_VAR(tstart, tfft1, taw, tfft2, tww, tr);
                                         //
                                         //    PROFILE_SAMPLE(tstart);

    //  Determine weighting filter spectrum W(exp(jw)) ---------------

    x[0] = ak[0];
    let mut coeff = gamma;
    for i in 1..order + 1 {
        x[i] = ak[i] * coeff;
        coeff *= gamma;
    }
    codec2_fftr(fftr_fwd_cfg, &x, &mut Ww);

    //    PROFILE_SAMPLE_AND_LOG(tfft2, taw, "        fft2");

    for i in 0..FFT_ENC / 2 {
        Ww[i].r = Ww[i].r * Ww[i].r + Ww[i].i * Ww[i].i;
    }

    //    PROFILE_SAMPLE_AND_LOG(tww, tfft2, "        Ww");

    //  Determined combined filter R = WA ---------------------------

    let mut max_Rw = 0.0;
    let mut min_Rw = 1E32;
    for i in 0..FFT_ENC / 2 {
        Rw[i] = (Ww[i].r * Pw[i]).sqrt();
        if Rw[i] > max_Rw {
            max_Rw = Rw[i];
        }
        if Rw[i] < min_Rw {
            min_Rw = Rw[i];
        }
    }

    //    PROFILE_SAMPLE_AND_LOG(tr, tww, "        R");

    //    #ifdef DUMP
    //    if (dump)
    //      dump_Rw(Rw);
    //    #endif

    //  create post filter mag spectrum and apply ------------------

    //  measure energy before post filtering

    let mut e_before = 1E-4;
    for i in 0..FFT_ENC / 2 {
        e_before += Pw[i];
    }
    //  apply post filter and measure energy

    //    #ifdef DUMP
    //    if (dump)
    //   dump_Pwb(Pw);
    //    #endif

    let mut e_after = 1E-4;
    for i in 0..FFT_ENC / 2 {
        let Pfw = Rw[i].powf(beta);
        Pw[i] *= Pfw * Pfw;
        e_after += Pw[i];
    }
    let mut gain = e_before / e_after;

    //  apply gain factor to normalise energy, and LPC Energy

    gain *= E;
    for i in 0..FFT_ENC / 2 {
        Pw[i] *= gain;
    }

    if bass_boost != 0 {
        //  add 3dB to first 1 kHz to account for LP effect of PF

        for i in 0..FFT_ENC / 8 {
            Pw[i] *= 1.4 * 1.4;
        }
    }

    //    PROFILE_SAMPLE_AND_LOG2(tr, "        filt");
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: interp_Wo()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22 May 2012

  Interpolates centre 10ms sample of Wo and L samples given two
  samples 20ms apart. Assumes voicing is available for centre
  (interpolated) frame.

\*---------------------------------------------------------------------------*/
pub fn interp_Wo(
    interp: &mut MODEL, //  interpolated model params
    prev: &MODEL,       //  previous frames model params
    next: &MODEL,       //  next frames model params
    Wo_min: f32,
) {
    interp_Wo2(interp, prev, next, 0.5, Wo_min);
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: interp_Wo2()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22 May 2012

  Weighted interpolation of two Wo samples.

\*---------------------------------------------------------------------------*/
fn interp_Wo2(
    interp: &mut MODEL, //  interpolated model params
    prev: &MODEL,       //  previous frames model params
    next: &MODEL,       //  next frames model params
    weight: f32,
    Wo_min: f32,
) {
    //  trap corner case where voicing est is probably wrong

    if interp.voiced != 0 && prev.voiced == 0 && next.voiced == 0 {
        interp.voiced = 0;
    }

    //  Wo depends on voicing of this and adjacent frames

    if interp.voiced != 0 {
        if prev.voiced != 0 && next.voiced != 0 {
            interp.Wo = (1.0 - weight) * prev.Wo + weight * next.Wo;
        }
        if prev.voiced == 0 && next.voiced != 0 {
            interp.Wo = next.Wo;
        }
        if prev.voiced != 0 && next.voiced == 0 {
            interp.Wo = prev.Wo;
        }
    } else {
        interp.Wo = Wo_min;
    }
    interp.L = (PI / interp.Wo as f64) as i32;
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: interpolate_lsp_ver2()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22 May 2012

  Weighted interpolation of LSPs.

\*---------------------------------------------------------------------------*/
pub fn interpolate_lsp_ver2(
    interp: &mut [f32],
    prev: &[f32],
    next: &[f32],
    weight: f32,
    order: usize,
) {
    for i in 0..order {
        interp[i] = (1.0 - weight) * prev[i] + weight * next[i];
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: apply_lpc_correction()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 22/8/2010

  Apply first harmonic LPC correction at decoder.  This helps improve
  low pitch males after LPC modelling, like hts1a and morig.

\*---------------------------------------------------------------------------*/
pub fn apply_lpc_correction(model: &mut MODEL) {
    if model.Wo < (PI as f32 * 150.0 / 4000.0) {
        model.A[1] *= 0.032;
    }
}

pub fn decode_lspds_scalar(lsp_: &mut [f32], indexes: &[usize], order: usize) {
    let mut lsp__hz = vec![0.0; order];
    let mut dlsp_ = vec![0.0; order];

    for i in 0..order {
        let k = lsp_cbd[i].k as usize;
        dlsp_[i] = lsp_cbd[i].cb[indexes[i] * k];

        if i != 0 {
            lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
        } else {
            lsp__hz[0] = dlsp_[0];
        }
        lsp_[i] = (PI as f32 / 4000.0) * lsp__hz[i];
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: lsp_to_lpc()
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 24/2/93

  This function converts LSP coefficients to LPC coefficients.  In the
  Speex code we worked out a way to simplify this significantly.

\*---------------------------------------------------------------------------*/
pub fn lsp_to_lpc(
    lsp: &[f32],
    ak: &mut [f32], // array of LPC coefficients
    order: usize,   // order of LPC coefficients
) {
    let mut freq = vec![0.0; order]; // array of LSP frequencies in radians
    let mut Wp = vec![0.0; (order * 4) + 2];

    //  convert from radians to the x=cos(w) domain
    for i in 0..order {
        freq[i] = lsp[i].cos();
    }

    let pwoff = 0; //Wp;
    let mut xin1 = 1.0;
    let mut xin2 = 1.0;

    //  reconstruct P(z) and Q(z) by cascading second order polynomials
    // in form 1 - 2xz(-1) +z(-2), where x is the LSP coefficient

    for j in 0..order + 1 {
        let mut n4off = 0;
        for i in 0..(order / 2) {
            let n1off = pwoff + (i * 4);
            let n2off = n1off + 1;
            let n3off = n2off + 1;
            n4off = n3off + 1;
            let xout1 = xin1 - 2.0 * (freq[2 * i]) * Wp[n1off] + Wp[n2off];
            let xout2 = xin2 - 2.0 * (freq[2 * i + 1]) * Wp[n3off] + Wp[n4off];
            Wp[n2off] = Wp[n1off];
            Wp[n4off] = Wp[n3off];
            Wp[n1off] = xin1;
            Wp[n3off] = xin2;
            xin1 = xout1;
            xin2 = xout2;
        }
        let xout1 = xin1 + Wp[n4off + 1];
        let xout2 = xin2 - Wp[n4off + 2];
        ak[j] = (xout1 + xout2) * 0.5;
        Wp[n4off + 1] = xin1;
        Wp[n4off + 2] = xin2;

        xin1 = 0.0;
        xin2 = 0.0;
    }
}

/*---------------------------------------------------------------------------*\

  FUNCTION....: synthesise
  AUTHOR......: David Rowe, conversion by Matt Weeks
  DATE CREATED: 20/2/95

  Synthesise a speech signal in the frequency domain from the
  sinusodal model parameters.  Uses overlap-add with a trapezoidal
  window to smoothly interpolate betwen frames.

\*---------------------------------------------------------------------------*/
pub fn synthesise(
    n_samp: usize,
    fftr_inv_cfg: &mut codec2_fftr_cfg,
    Sn_: &mut [f32], //  time domain synthesised signal
    model: &MODEL,   //  ptr to model parameters for this frame
    Pn: &[f32],      //  time domain Parzen window
    shift: bool,     //  flag used to handle transition frames
) {
    let mut Sw_ = [COMP::new(); FFT_DEC / 2 + 1]; //  DFT of synthesised signal
    let mut sw_ = [0.0; FFT_DEC]; //  synthesised signal

    if shift {
        //  Update memories
        for i in 0..n_samp - 1 {
            Sn_[i] = Sn_[i + n_samp];
        }
        Sn_[n_samp - 1] = 0.0;
    }

    for i in 0..FFT_DEC / 2 + 1 {
        Sw_[i].r = 0.0;
        Sw_[i].i = 0.0;
    }

    //  Now set up frequency domain synthesised speech
    for l in 1..model.L as usize + 1 {
        let mut b = (l as f32 * model.Wo * FFT_DEC as f32 / TWO_PI + 0.5) as usize;
        if b > ((FFT_DEC / 2) - 1) {
            b = (FFT_DEC / 2) - 1;
        }
        Sw_[b].r = model.A[l] * model.phi[l].cos();
        Sw_[b].i = model.A[l] * model.phi[l].sin();
    }

    //  Perform inverse DFT
    codec2_fftri(fftr_inv_cfg, &Sw_, &mut sw_);

    //  Overlap add to previous samples

    //    #ifdef USE_KISS_FFT
    //    #define    FFTI_FACTOR ((float)1.0)
    //    #else
    //    #define    FFTI_FACTOR ((float32_t)FFT_DEC)
    //    #endif

    for i in 0..n_samp - 1 {
        Sn_[i] += sw_[FFT_DEC - n_samp + 1 + i] * Pn[i]; // * FFTI_FACTOR;
    }

    if shift {
        let mut j = 0;
        for i in n_samp - 1..2 * n_samp {
            Sn_[i] = sw_[j] * Pn[i]; // * FFTI_FACTOR;
            j += 1;
        }
    } else {
        let mut j = 0;
        for i in n_samp - 1..2 * n_samp {
            Sn_[i] += sw_[j] * Pn[i]; // * FFTI_FACTOR;
            j += 1;
        }
    }
}

const BG_THRESH: f32 = 40.0; //  only consider low levels signals for bg_est
const BG_BETA: f32 = 0.1; //  averaging filter constant
const BG_MARGIN: f32 = 6.0; //  harmonics this far above BG noise are
                            //   randomised.  Helped make bg noise less
                            //   spikey (impulsive) for mmt1, but speech was
                            //               perhaps a little rougher.

/*---------------------------------------------------------------------------*\

  postfilter()

  The post filter is designed to help with speech corrupted by
  background noise.  The zero phase model tends to make speech with
  background noise sound "clicky".  With high levels of background
  noise the low level inter-formant parts of the spectrum will contain
  noise rather than speech harmonics, so modelling them as voiced
  (i.e. a continuous, non-random phase track) is inaccurate.

  Some codecs (like MBE) have a mixed voicing model that breaks the
  spectrum into voiced and unvoiced regions.  Several bits/frame
  (5-12) are required to transmit the frequency selective voicing
  information.  Mixed excitation also requires accurate voicing
  estimation (parameter estimators always break occasionally under
  exceptional conditions).

  In our case we use a post filter approach which requires no
  additional bits to be transmitted.  The decoder measures the average
  level of the background noise during unvoiced frames.  If a harmonic
  is less than this level it is made unvoiced by randomising it's
  phases.

  This idea is rather experimental.  Some potential problems that may
  happen:

  1/ If someone says "aaaaaaaahhhhhhhhh" will background estimator track
     up to speech level?  This would be a bad thing.

  2/ If background noise suddenly dissapears from the source speech does
     estimate drop quickly?  What is noise suddenly re-appears?

  3/ Background noise with a non-flat sepctrum.  Current algorithm just
     comsiders scpetrum as a whole, but this could be broken up into
     bands, each with their own estimator.

  4/ Males and females with the same level of background noise.  Check
     performance the same.  Changing Wo affects width of each band, may
     affect bg energy estimates.

  5/ Not sure what happens during long periods of voiced speech
     e.g. "sshhhhhhh"

\*---------------------------------------------------------------------------*/
pub fn postfilter(model: &mut MODEL, bg_est: &mut f32) {
    //  determine average energy across spectrum

    let mut e = 1E-12;
    for m in 1..model.L as usize + 1 {
        e += model.A[m] * model.A[m];
    }
    //  assert(e > 0.0);
    e = 10.0 * (e / model.L as f32).log10();

    //  If beneath threhold, update bg estimate.  The idea
    // of the threshold is to prevent updating during high level
    // speech.

    if e < BG_THRESH && model.voiced != 0 {
        *bg_est = *bg_est * (1.0 - BG_BETA) + e * BG_BETA;
    }
    //  now mess with phases during voiced frames to make any harmonics
    //  less then our background estimate unvoiced.

    //  let mut uv = 0;
    let thresh = 10.0_f32.powf((*bg_est + BG_MARGIN) / 20.0);
    if model.voiced != 0 {
        for m in 1..model.L as usize + 1 {
            if model.A[m] < thresh {
                model.phi[m] = (TWO_PI / CODEC2_RAND_MAX) * codec2_rand();
                //uv++;
            }
        }
    }
    //#ifdef DUMP
    //  dump_bg(e, *bg_est, 100.0*uv/model.L);
    //#endif
}
