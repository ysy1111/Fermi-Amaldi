This file records all the annotations/changes/modifications in ~/nwchem-7.0.2-release source code

#Annotation style:
    ---start with: c=ysy/begin=c
            code here
    ---end with: c=ysy/end=c
    one template:
    
    c=ysy/begin=c
        OPEN(10, file='info.txt', status='old',access='append')
        WRITE(10,*) '.F:'
        WRITE(10,*)
        CLOSE(10)
    c=ysy/end=c
#record style:
    ---path/file_name(line number)
    ---usage of the code
    ---anything important

1.dft_scf.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/scf_dft/dft_scf.F
---
variable identification: use_nwxc = True
==============================================================
2.nwxc_query.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwxc/nwxc_query.F
---
xc_related, has subroutine nwxc_query and nwxc_getvals
===============================================================
3.nwxcP.fh
path: /home/shuyy/nwchem-7.0.2-release/src/nwxc/nwxcP.fh
---
xc related: identify the functionals' names
================================================================
4.xc_getv.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/xc/xc_getv.F
---
Calculate the exchange-correlation energy and Fock matrix
Call xc_getvxc.F to compute the dft xc
---
call fock_2e to build fock matrix? need g_dens, g_jk to calculate coulomb energy
---
The energy part can be directly modified here.
---
Check the potential part, how to implement the FA potential
---
===================================================================
5.xc_getvxc.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/xc/xc_getvxc.F
---
call subroutine grid_quadv0-->further call grid_quadv0_gen
====================================================================
6.grid_quadv0.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/grid/grid_quadv0.F
---
subroutine grid_quadv0_gen ### important ###
====================================================================
7.grid_quadv0a.F,grid_quadvw.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/grid/grid_quadv0a.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/grid/grid_quadvw.F
---
after these two subroutine the energy has been calculated.<--check how
---
grid_quadvw.F -> call grid_loop
---
grid_quadv0a.F -> call grid_quadv0b --> call xc_eval_fnl
====================================================================
8.xc_eval_fnl.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/xc/xc_eval_fnl.F
---
functional evaluation routine
====================================================================
9.
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/xc/xc_dirac
---
calculation dirac exchange energy
====================================================================
10.fock_2e.F
path:/home/shuyy/nwchem-7.0.2-release/src/ddscf/fock_2e.F
---
called by xc_getv.F to build fock matrix?
=======================================================
11.ao_fock_2e->ao_replicated
path:/home/shuyy/nwchem-7.0.2-release/src/ddscf/ao_fock_2e.F
path:/home/shuyy/nwchem-7.0.2-release/src/ddscf/ao_replicated.F
---
subroutine ao_fock_2e is called by fock_2e
subroutine ao_replicated is called by ao_fock_2e
inside ao_replicated.F, subroutine ao_replicated calls subroutine fock_2e_rep for the first two rounds then fock_2e_rep_from_file from fock_2e_file.F three rounds.
---
==============================================
12.xc_tabcd.F
path:/home/shuyy/nwchem-7.0.2-release/src/nwdft/xc/xc_tabcd.F
==============================================
13.fock_2e_file.F: -> fock_2e_rep_from_file
path:/home/shuyy/nwchem-7.0.2-release/src/ddscf
---
call int2e_file_rep_fock? to build update fock matrices?



====================================================
=====================*WORK FLOW*====================
===================*Unrestricted*===================
====================================================


src/nwdft/scf_dft/dft_scf.F(main scalar DFT driver)
|
--> dft_fockbld ==> build the fock matrix
           |
           -->xc_getv ==> calculate the HF exchange/Coulomb & DFT xc
                   |
                   --> fock_2e ==> calculate the coulomb parts
                            |
                            --> ao_fock_2e 
 |<-------------------------------------|
 |
 --> ao_replicated------------>|
              |                \
              --> fock_2e_rep   |-->fock_2e_rep_from_file(fock_2e_file.F)
                  (first)          (later)
  |<------------------|            |     
  |                                |
  -->fock_rep_txs(ao_replicated.F) -->int2e_file_rep_fock(int2e_file.F)
               |                       |
*fock_2e_lab.F*-->fock_2e_rep_label    -->fock_2e_rep_mod_label ------->|
 |               |                                                      |
 |               -->fock_2e_rep_1_label(before build fock)**            |
 |               |                                                      |
 |               -->fock_2e_rep_4_label(during building)**              |
 |                                                                      |
 |               |<-----------------------------------------------------|
 |               |
                 --> fock_2e_rep_mod_4_label****


====================================================
=====================*WORK FLOW*====================
====================*Restricted*====================
====================================================


Version 1:

    jfac --> jfac * (1-a2/N)

    potential:
    j_a --> (1-a2/N)j_a
    j_b --> (1-a2/N)j_b
    j = j_a + j_b --> (1-a2/N)j
    f_a,f_b both contain j --> (1-a2/N)j

    Energy:
    J_a --> (d_a,j) = (d_a,j_a) + (d_a,j_b)
    J_b --> (d_b,j) = (d_b,j_a) + (d_b,j_b)
    J = J_a + J_b =(d,j)=(1-a2/N)J

Version 2:

    k_a -->  a2/N_a * j_a
    k_b -->  a2/N_b * j_b

    f_a --> k_a
    f_b --> k_b

    K_a --> -(d_a,k_a) =  -a2/N_a * J_a
    K_b --> -(d_b,k_b) =  -a2/N_b * J_b

    if N_a == N_b == N/2:
    k_a == k_b == -a2/(N/2) * j/2 = -a2/N * j
    K_a = -(d_a,k_a) = -(d/2,k_a) = -a2/(2N) * J
    K_b = -(d_b,k_b) = -(d/2,k_b) = -a2/(2N) * J
    K = K_a+K_b = -a2/N * J