use crate::*;

/* codebook/dlsp1.txt */
const codes0: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp2.txt */
const codes1: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp3.txt */
const codes2: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp4.txt */
const codes3: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
    550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
    1200.0, 1250.0, 1300.0, 1350.0, 1400.0,
];
/* codebook/dlsp5.txt */
const codes4: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
    550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
    1200.0, 1250.0, 1300.0, 1350.0, 1400.0,
];
/* codebook/dlsp6.txt */
const codes5: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0,
    550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
    1200.0, 1250.0, 1300.0, 1350.0, 1400.0,
];
/* codebook/dlsp7.txt */
const codes6: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp8.txt */
const codes7: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp9.txt */
const codes8: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];
/* codebook/dlsp10.txt */
const codes9: [f32; 32] = [
    25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0,
    375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 725.0, 750.0, 775.0, 800.0,
];

const codesnone: [f32; 0] = [];

pub const lsp_cbd: [lsp_codebook; 11] = [
    /* codebook/dlsp1.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes0,
    },
    /* codebook/dlsp2.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes1,
    },
    /* codebook/dlsp3.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes2,
    },
    /* codebook/dlsp4.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes3,
    },
    /* codebook/dlsp5.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes4,
    },
    /* codebook/dlsp6.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes5,
    },
    /* codebook/dlsp7.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes6,
    },
    /* codebook/dlsp8.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes7,
    },
    /* codebook/dlsp9.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes8,
    },
    /* codebook/dlsp10.txt */
    lsp_codebook {
        k: 1,
        log2m: 5,
        m: 32,
        cb: &codes9,
    },
    lsp_codebook {
        k: 0,
        log2m: 0,
        m: 0,
        cb: &codesnone,
    },
];
