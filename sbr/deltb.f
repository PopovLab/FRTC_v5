       SUBROUTINE DELTB(D0, U0, I0, R2)
       INCLUDE 'for/parameter.inc'
       INCLUDE 'for/const.inc'
       INCLUDE 'for/status.inc'

       REAL*8 B, U0, I0, TINTG, D0, R2, BTR, LNR

       B = 0.103
       TINTG = 0.00136
       CV1 = -B*TINTG*(WPFL-U0) / 0.2 / GP / RTOR / IPL + D0*I0 / IPL

       YB = 0.
       do 1 J=1, NA
       YB = YB + ELON(J)*AMETR(J)**2 * 
     & (NE(J)*TE(J) + NI(J)*TI(J) - NE(J+1)*TE(J+1) - NI(J+1)*TI(J+1))
1      CONTINUE

       BTR = 0.00064*GP2*YB*
     & (RTOR/(G22G(NA)*IPOL(NA)*BTOR*NA*HRO*MU(NA)))**2
       BTR = BTR*R2*R2/AB/AB

       SINT = 0
       do 2 J=1,NA
       J1 = J + J - 1
2      SINT = SINT + (J1 * MU(J)**2 * (2*J*J-J1)) * IPOL(J) * G22G(J)/J
       LNR = 0.5 * SINT * HRO * 
     & ((NA+1)/(NA*NA*MU(NA+1)*IPOL(NA+1)*G22G(NA+1)))**2
       LNR = LNR + 2*LOG(R2/AB)

       CV2 = CV1 + B*B*LOG(B/R2)/2/RTOR
       CV2 = CV2 + (B*B - R2*R2) * (BTR + (LNR-1)/2)/2/RTOR

       END