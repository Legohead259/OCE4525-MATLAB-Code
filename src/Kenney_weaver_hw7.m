%% Owner: Clay Kenney
    %Date: 11/22/2021
    %Student ID #: 903577008

    %File Name: Kenney_weaver_hw7.m
    
    clear all
    close all
    
%% Inputs
    h = 6;
    bottom_slope = .02;
    H_0 = 5.8;
    T = 20;
    theta_0 = 0;
    rho = 1025;
    g = 9.91;
    tolerance = .00001;

%% Calculations
    %Part A Section 1
        rho_s = 2650;
        rho_w = rho_s - rho;
        delta = rho_s / rho_w;
        a = .013;
        b = 22;
        Hs = H_0;
        cota = 2;
        Kd1 = 2.0;
        M50_1 = (rho_s * Hs ^ 3) / (Kd1 * (delta ^ 3) * cota);
        Dn50_1 = (M50_1 / rho_s)^(1/3);
        n1 = 2;
        kdell = 1;
        r1 = n1 * kdell * Dn50_1;
        P1 = 37;
        gamma_1 = rho_s * g;
        NaN1 = n1 * kdell * (rho_s/gamma_1)^(2/3) * (1 - P1/100);
        [L, L_0, k, kh, C, Cg, k_0, kh_0, C_0, Cg_0, Kr, Ks, theta_A, E, F, R, sigma, H, H_0_prime] = kenney_wave_reader(T, h, H_0, theta_0, g, rho, tolerance);
        L0 = L_0;
        tanb = .03;
        delta = (rho_s / rho_w) - 1;
        R1 = Hs / (delta*Dn50_1);
        Ns1 = R1 * (Hs / L)^(-1/3);
        etam1 = Hs * tanb / sqrt(Hs / L0);
        
        Hs25 = H_0 / 1.37;
        Kd2 = 31;
        cota2 = 2;
        M50_2 = (rho_s * Hs25 ^ 3) / (Kd2 * (delta ^ 3) * cota2);    
        Dn50_2 = (M50_2 / rho_s)^(1/3);
        r2 = n1 * kdell * Dn50_2;
        P2 = P1;
        W50 = M50_2 * g;
        NaN2 = n1 * kdell * (gamma_1/W50)^(2/3) * (1 - P2/100);
        etam2 = Hs * tanb / sqrt(Hs25 / L0);

    %% Part A Section 2        
        rho_s = 2400;
        delta3 = (rho_s / rho)-1;        
        Hs3 = H_0;
        cota = 2;
        Kd3 = 7;
        M50_3 = (rho_s * Hs3 ^ 3) / (Kd3 * (delta3 ^ 3) * cota);
        Dn50_3 = (M50_3 / rho_s)^(1/3);
        kdel3 = 1.04;
        r2 = n1 * kdel3 * Dn50_3;
        W503 = M50_3 * g;
        gamma_2 = rho_s * g;
        P3 = 50;
        NaN2 = n1 * kdel3 * (gamma_2/W503)^(2/3) * (1 - P3/100);
        etam3 = Hs3 * tanb / sqrt(Hs3 / L0);
        
        Hs3_25 = Hs / 1.32;
        Kd3 = 7;
        cota = 2;
        M50_4 = (rho_s * Hs3_25 ^ 3) / (Kd3 * (delta3 ^ 3) * cota);    
        Dn50_4 = (M50_4 / rho_s)^(1/3);
        r4 = n1 * kdel3 * Dn50_4;
        W50 = M50_4 * g;
        NaN4 = n1 * kdel3 * (gamma_2/W50)^(2/3) * (1 - P3/100);
        etam4 = Hs * tanb / sqrt(Hs3_25 / L0);        
        

    %% Part B Section 1
        Nz = 7500;
        P5 = 0.4;
        rho_s = 2650;
        rho_w = rho_s - rho;
        delta = (rho_s / rho_w)-1
        Hs = H_0;
        cota = 2;
        Kd5 = 24;
        M50_5 = (rho_s * Hs ^ 3) / (Kd5 * (delta ^ 3) * cota);
        Dn50_5 = (M50_5 / rho_s)^(1/3);
        r5 = n1 * kdell * Dn50_5;
        sm5 = Hs / L0;
        S = 2;
        Ae5 = S * Dn50_5^2;

        delta= (rho_s / rho)-1
        Ns = Hs / (delta * Dn50_5)
        NaN5 = Ns / Ae5        

        
        
        
        
        
        