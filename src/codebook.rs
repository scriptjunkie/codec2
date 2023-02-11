use crate::*;

/* THIS IS A GENERATED FILE. Edit generate_codebook.c and its input */

/*
 * This intermediary file and the files that used to create it are under
 * The LGPL. See the file COPYING.
 */

/* codebook/lsp1.txt */
const codes0: [f32; 16] = [
    225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0,
    550.0, 575.0, 600.0,
];
/* codebook/lsp2.txt */
const codes1: [f32; 16] = [
    325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0,
    650.0, 675.0, 700.0,
];
/* codebook/lsp3.txt */
const codes2: [f32; 16] = [
    500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0,
    1150.0, 1200.0, 1250.0,
];
/* codebook/lsp4.txt */
const codes3: [f32; 16] = [
    700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0,
    1900.0, 2000.0, 2100.0, 2200.0,
];
/* codebook/lsp5.txt */
const codes4: [f32; 16] = [
    950.0, 1050.0, 1150.0, 1250.0, 1350.0, 1450.0, 1550.0, 1650.0, 1750.0, 1850.0, 1950.0, 2050.0,
    2150.0, 2250.0, 2350.0, 2450.0,
];
/* codebook/lsp6.txt */
const codes5: [f32; 16] = [
    1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0,
    2300.0, 2400.0, 2500.0, 2600.0,
];
/* codebook/lsp7.txt */
const codes6: [f32; 16] = [
    1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0,
    2700.0, 2800.0, 2900.0, 3000.0,
];
/* codebook/lsp8.txt */
const codes7: [f32; 8] = [
    2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0,
];
/* codebook/lsp9.txt */
const codes8: [f32; 8] = [
    2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0,
];
/* codebook/lsp10.txt */
const codes9: [f32; 4] = [2900.0, 3100.0, 3300.0, 3500.0];

pub const lsp_cb: [lsp_codebook; 10] = [
    /* codebook/lsp1.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes0,
    },
    /* codebook/lsp2.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes1,
    },
    /* codebook/lsp3.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes2,
    },
    /* codebook/lsp4.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes3,
    },
    /* codebook/lsp5.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes4,
    },
    /* codebook/lsp6.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes5,
    },
    /* codebook/lsp7.txt */
    lsp_codebook {
        k: 1,
        log2m: 4,
        m: 16,
        cb: &codes6,
    },
    /* codebook/lsp8.txt */
    lsp_codebook {
        k: 1,
        log2m: 3,
        m: 8,
        cb: &codes7,
    },
    /* codebook/lsp9.txt */
    lsp_codebook {
        k: 1,
        log2m: 3,
        m: 8,
        cb: &codes8,
    },
    /* codebook/lsp10.txt */
    lsp_codebook {
        k: 1,
        log2m: 2,
        m: 4,
        cb: &codes9,
    },
];
