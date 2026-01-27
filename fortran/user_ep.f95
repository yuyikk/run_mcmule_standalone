!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE USER
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use MCMULE
    implicit none

    integer, parameter :: nrq = 1
    integer, parameter :: nrbins = 750
    real(kind = prec), parameter :: deg = 180. / pi
    integer :: process
    real(kind = prec) :: e_beam

    real(kind = prec), parameter :: min_val(nrq) = (/0.5/)
    real(kind = prec), parameter :: max_val(nrq) = (/8.0/)

    integer :: userdim = 0


    integer, parameter :: namesLen=6
    integer, parameter :: filenamesuffixLen=10
    integer :: nq=nrq
    integer :: nbins=nrbins



    integer :: bin_kind = 0;

    contains


    SUBROUTINE FIX_MU

        musq = me**2

    
    END SUBROUTINE FIX_MU

    SUBROUTINE INITUSER
        print*, "PRad with electrons"
        print*, " * dipole with Lambda^2 = 0.88 GeV^2"
        print*, " * electron energy beam E_e = 1100 MeV"
        print*, " * angular window: 0.5 < theta_e < 8. deg" 

        print*, " * Please enter beam energy (MeV): "
        read*, e_beam
        print*, " * running ep-scattering"
        print*, " * additional cut on photons with E_y > 20 MeV if th(el_f,y)>6 mrad"
        if(which_piece(1:5) == "mp2mp") then
            call initflavour("e-p-", Mel**2 + Mproton**2 + 2 * Mproton * e_beam)
        print*, " * running QED-FF with e^-"
        endif
    END SUBROUTINE

    FUNCTION QUANT(q1, q2, q3, q4, q5, q6, q7)

        real (kind=prec), intent(in) :: q1(4),q2(4),q3(4),q4(4), q5(4),q6(4),q7(4)
        real (kind=prec) :: q1Rest(4),q2Rest(4),q3Rest(4),q4Rest(4),q5Rest(4),q6Rest(4)
        real (kind=prec) :: quant(nrq)
        real (kind=prec) :: th_l

        pol1 = (/ 0._prec, 0._prec, 0._prec, 0._prec /)

        pass_cut = .true.
        call FIX_MU

        q1Rest = boost_rf(q2,q1)  ! incomping lepton
        q2Rest = boost_rf(q2,q2)  ! lepton/proton at rest
        q3Rest = boost_rf(q2,q3)  ! outgoing lepton
        q4Rest = boost_rf(q2,q4)  ! recoiling lepton/proton
        q5Rest = boost_rf(q2,q5)  ! outgoing photon (if present)
        q6Rest = boost_rf(q2,q6)  ! outgoing photon (if present)

        th_l = acos(cos_th(q1Rest, q3Rest))

        names(1) = "th_l"
        quant(1) = th_l
    END FUNCTION QUANT