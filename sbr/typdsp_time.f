C Modules >>> TYPDSP, PCOVA, UF1DWA, UF2DWA, OPENAP, OPENRD, OPENWT <<<
C >>> SETFRM, AFRAME, DNSTR, UPSTR,  TIMEDT, PUTXY,  NEGA,   SETFNA <<< 
C >>> SKIPM  <<< 
C======================================================================|
	subroutine	TYPDSP_TIME(NCH,YN,ITIMES,TTOUT,TOUT)
C NCH= 5 - terminal, 0 - file (old format), 1 - file (new format)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	logical, save :: first_time=.TRUE.
	integer	ITIMES
	double precision	YN,YQ,TTOUT(ITIMES),TOUT(ITIMES,NRW)
	integer NCH,NP1,ITBE,ITEND,ITEN,IERR,MODEX,NLINSC,NVAR,NNJ,JN0
	integer JBE,JEND,J,JEN,JJ,J1,JLR,KILLBL,length,WarningColor
	character*6 CH6,STRMN*40,CONN(6)
	character(132) FNAME, STRI, tmp
	logical	EXI

	FNAME = "dat/time_series.dat"
	JEND = KILLBL(FNAME,132)
	print '(">>>  time series are written in: ", A)', FNAME(1:JEND)

		open(20, file=FNAME,position="append")
	
		if (first_time) then
			write(20,'(100A22)'), "Time", (NAMET(J),J=1,NTOUT)
			first_time = .FALSE.
		  end if
		write(20,'(100f22.14)') TIME, (TOUT(LTOUT,j), j=1,NTOUT )
    	close(20)
	END

