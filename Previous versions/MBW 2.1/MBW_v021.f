C -----------------------------------------------------------C
C Subroutine for Abaqus/Explicit for isotropic elasticity    C
C and isotropic plasticity with pressure and lode            C
C dependence, according to plasticity model by BAI, Y. and   C
C WIERZBICKI, T. (published in "A new model of metal         C
C plasticity and fracture with pressure and Lode             C
C dependence", International Journal of Plasticity 24 (2008) C
C pages 1071-1096).                                          C
C															 C
C - Damage initiation according to 3D fracture locus         C
C   (equivalent plastic strain to failure initiation being   C
C   function of both hydrostatic pressure and lode angle)    C
C - Damage evolution based on energy dissipated during       C
C   damage process                                           C
C															 C
C Originally based on the von Mises' J2 plasticity theory.   C
C Elastic predictor, radial corrector algorithm used with    C
C implicit backward Euler scheme for integration of flow     C
C rule and damage evolution equation.                        C
C                                                            C
C Not suitable for solving problem settings including plane  C
C stresses. Not suitable for 2D models.                      C
C -----------------------------------------------------------C
C Solution dependent variables (SDVs):						 C
C															 C
C SDV1: 	Equivalent plastic strain						 C
C SDV2: 	Damage variable									 C
C SDV3: 	Yield stress at damage initiation				 C
C SDV4: 	Flag (=0 element not damaged,					 C
C			      =1 element experienced brittle damage,	 C
C			      =2 element experienced ductile damage)	 C
C SDV5-10:	Total strain tensor								 C
C SDV11-16:	Plastic strain tensor							 C
C SDV17:	Equivalent plastic strain increment				 C
C SDV18:	Temperature softening correction function		 C
C SDV19:	Mutliplicative strain rate hardening correction	 C
C			function		 								 C
C SDV20:	Equivalent plastic strain rate					 C
C SDV21:	Equivalent plastic strain rate (computed over	 C
C           a user-defined no. of increments)				 C
C SDV22:	Increment counter (for eq. plas. str. rate       C
C           computation)									 C
C SDV23:	Sum of equivalent plastic strain increments      C
C           before the computation of eq. plas. str. rate    C
C SDV24:	Sum of time increments before computation of eq. C
C           plas. str. rate									 C
C SDV25:	Additive strain rate hardening correction		 C
C			function										 C
C SDV26:	Maximum principal stress						 C
C SDV27:	Equivalent plastic strain (averaged over a		 C
C			user-defined no. of increments)					 C
C SDV28:	Element deletion flag							 C
C SDV29:	Equivalent plastic strain to failure initiation  C
C -----------------------------------------------------------C
C Contact info:                                              C
C															 C
C Mohamed Sharaf                                             C
C RWTH Aachen                                                C
C Institut fuer Eisenhuettenkunde                            C
C Intzestrasse 1                                             C
C 52072 Aachen                                               C
C Germany                                                    C
C mohamed.sharaf@iehk.rwth-aachen.de                         C
C -----------------------------------------------------------C
C Aachen, 19 March 2012                                      C
C -----------------------------------------------------------C
        subroutine vumat(
C Read only
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     1  stepTime, totalTime, dt, cmname, coordMp, charLength,
     1  props, density, strainInc, relSpinInc,
     1  tempOld, stretchOld, defgradOld, fieldOld,
     1  stressOld, stateOld, enerInternOld, enerInelasOld,
     1  tempNew, stretchNew, defgradNew, fieldNew,
C Write only
     1  stressNew, stateNew, enerInternNew, enerInelasNew)
C
        include 'vaba_param.inc'
C
        dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     1  relSpinInc(*), tempOld(nblock),
     1  stretchOld(*), defgradOld(*),
     1  fieldOld(*), stressOld(nblock,ndir+nshr),
     1  stateOld(nblock,nstatev), enerInternOld(nblock),
     1  enerInelasOld(nblock), tempNew(nblock),
     1  stretchNew(*), defgradNew(*), fieldNew(*),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
        character*80 cmname
C
C Defining numerical constants
        parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     1  threeHalves = 1.5, thousand = 1000. )
        data newton,tol/10,1.D-6/
        data charLengthh,damold/0.D0,0.D0/
        data dam/0.D0/
C
C
C Defining material properties constants
        e      = props(1)
        xnu    = props(2)
        ceta   = props(3)
        eta0   = props(4)
        cthetas= props(5)
        cthetat= props(6)
        cthetac= props(7)
        om     = props(8)
C       Ductile damage initiation locus function parameters
        d1     = props(9)
        d2     = props(10)
        d3     = props(11)
        d4     = props(12)
        d5     = props(13)
        d6     = props(14)
C       Fracture energy dissipation (energy dissipated per element unit area after
C       ductile damage dissipation) [Joules/mm2]
        gf     = props(15)
C       Maximum value for ductile damage variable
        Dcrit   = props(16)
C       5 parameters for temperature softening correction function
        beta   = props(17)
        alpha  = props(18)
        eta2   = props(19)
        cp     = props(20)
        t0     = props(21)
C       2 parameters for strain rate hardening correction function
        p2     = props(22)
        d     = props(23)
C       Number of increments to update plastic strain rate
        strrInc= props(24)
C       Brittle damage initiation stress value
        sigdmg = props(25)
C
        nvalue = (nprops-32)/2
C
C Checking for valid entries
        if (om.lt.zero) then
          write(6,5)
 5        format(//,30X,'***ERROR - m MUST BE A NON-NEGATIVE INTEGER')
        endif
C
C Defining constants
        twomu  = e / ( one + xnu )
        thremu = threeHalves * twomu
        sixmu  = three * twomu
        alamda = twomu * ( e - twomu ) / ( sixmu - two * e )
        akappa = e * twomu * half / (three * (threeHalves * twomu - e 
     1           ))
        term   = one / ( twomu * ( one + hard/thremu ) )
        con1   = sqrt( twoThirds )
        pi     = 4. * atan(one)
        con2   = (sqrt(three)/two) / (one - sqrt(three)/two)
C
C Iteration starts here
      do 100 i = 1,nblock
C
C If brittle damage has previously occured, stress is zeroed
        if (stateOld(i,4).gt.half .and. stateOld(i,4).lt.(one+half))
     1    then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4) = one
          goto 100
        endif
C Updating total strain (states 5-10)
        stateNew(i,5)  = stateOld(i,5)  + strainInc(i,1)
        stateNew(i,6)  = stateOld(i,6)  + strainInc(i,2)
        stateNew(i,7)  = stateOld(i,7)  + strainInc(i,3)
        stateNew(i,8)  = stateOld(i,8)  + strainInc(i,4)
        stateNew(i,9)  = stateOld(i,9)  + strainInc(i,5)
        stateNew(i,10) = stateOld(i,10) + strainInc(i,6)
C
C       Trace of total strain tensor
        epsilontrace = stateNew(i,5) + stateNew(i,6) + stateNew(i,7)
C
C       Trace of strain increment tensor
        trace = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
C
C Calculating softening correction term due to temperature rise
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
C          facT OFF
C          facT = beta / t0         **alpha
          facT = one
          stateNew(i,18)= facT
        else
C          facT OFF
C          facT = beta / tempOld(i) **alpha
          facT = one
          stateNew(i,18)= facT
        endif   
C
C Calculating hardening correction term due to straining rate
C        facStrrate OFF
C        facStrrate = p2 * log(stateOld(i,21)) + d
        facStrrate = one
C
C Preventing negative values for straining rate
        if (facStrrate.lt. one) facStrrate = one
        stateNew(i,19)= facStrrate
C
C Restoring previous softening correction term due to ductile damage
C        DAMAGE ON
        facDctlDmgDev = one - stateOld(i,2)
C        facDctlDmgDev = one
        facDctlDmgVol = facDctlDmgDev
        if (stressOld(i,1)+stressOld(i,2)+stressOld(i,3).lt.zero) then
          facDctlDmgVol = one
        endif
C
C Trial stress
        sig1 = stressOld(i,1) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,1)
        sig2 = stressOld(i,2) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,2)
        sig3 = stressOld(i,3) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,3)
        sig4 = stressOld(i,4) + facDctlDmgDev*twomu*strainInc(i,4)
        sig5 = stressOld(i,5) + facDctlDmgDev*twomu*strainInc(i,5)
        sig6 = stressOld(i,6) + facDctlDmgDev*twomu*strainInc(i,6)
C
C Following equation numbering of publication mentioned in code title text
C Equ. (4) - Calculating the deviatoric part of trial stress
        sigmean = third * ( sig1 + sig2 + sig3 )
        ds1 = sig1 - sigmean
        ds2 = sig2 - sigmean
        ds3 = sig3 - sigmean
C
C Calculating the magnitude of the deviatoric trial stress tensor
        dsmag = sqrt( ds1**2 + ds2**2 + ds3**2 + two*sig4**2 + two
     1  *sig5**2
     1  + two*sig6**2 )
C
C Preventing a divide by zero when dsmag is zero. When dsmag is zero, computation is still in the elastic zone
        if (dsmag.lt.tol .and. dsmag.ge.zero) dsmag = tol
        if (dsmag.gt.(zero-tol) .and. dsmag.le.zero) dsmag = zero-tol
C
C Following equation numbering of publication mentioned in code title text
C Eq. (1) - Calculating the 1st invariant of the stress tensor
        p   = zero - sigmean
C
C Eq. (2) - Calculating the 2nd invariant of the stress tensor
        q   = dsmag / con1
C
C Eq. (3) - Calculating the 3rd invariant of the stress tensor
        r   = ( three/twoThirds * (ds1*(ds1**2 + sig4**2 + sig6**2)
     1  + 2.
     1  *sig4*(ds1*sig4 + ds2*sig4 + sig6*sig5) + 2.*sig6*(ds1*sig6
     1  + sig4
     1  *sig5 + ds3*sig6) + ds2*(sig4**2 + ds2**2 + sig5**2) + 2.*sig5
     1  *(sig4*sig6 + ds2*sig5 + ds3*sig5) + ds3*(sig6**2 + sig5**2
     1  + ds3**2)) )**third
C
C Eq. (5) - Calculating the dimensionless hydrostatic pressure eta
        eta = sigmean / q
C
C Eq. (6) - Calculating the Lode angle theta
        cosine=(r/q)**3
C       Assuring that -1<cosine<1
        if (cosine.gt.one) cosine=one
        if (cosine.lt.(zero-one)) cosine=zero-one
        theta = third*acos(cosine)
C
C Eq. (7) - Calculating the normalized Lode angle thetabar
        thetabar = one - 6.*theta/pi
C
C Eq. (12) - Calculating the Lode angle dependant parameter gamma
        gamma = con2 * (one/cos(theta-pi/6.) - one)
C
C Eq. (13) - Determining cthetaax, either tension case or compression case
        if( thetabar .ge. zero ) then
          cthetaax = cthetat
        else
          cthetaax = cthetac
        endif
C
C Fetching yield stress from flow curve
        call ahard(sigmayield,hard,stateOld(i,1),props(33),nvalue)
C
C Eq. (11) - Calculating radius of elastic zone
        fita   = one - ceta*(eta-eta0)
        ftheta = cthetas + ((cthetaax-cthetas)* (gamma- ((gamma**(om
     1  +one))/(om
     1  +one)) ) ) 
        radius = sigmayield * fita * ftheta * facT * facStrrate
     1  * facDctlDmgVol
C
C Eq. (15) - Checking if yielding. The yielding factor facyld is set to be zero
C in the elastic zone, one if yielding
        facyld = zero
        if( q - radius .ge. zero ) facyld = one
C
C Eq. (20) - Calculating the tensor of the normal direction of plastic flow
C       Avoiding a divide by zero by avoiding that theta = zero
        if (theta.eq.zero) theta = theta + tol
C       Derivative of triaxiality correction function w.r.t. stress
        dfitads  = ceta * three*eta/(two*(q**2))
        dfitads1 = dfitads * ds1
        dfitads2 = dfitads * ds2
        dfitads3 = dfitads * ds3
        dfitads4 = dfitads * sig4
        dfitads5 = dfitads * sig5
        dfitads6 = dfitads * sig6
C       Derivative of Lode angle correction function w.r.t. stress
        dfthetads  = (cthetaax-cthetas) * (one -(gamma**om)) * (three
     1  *sqrt(three)/(two-sqrt(three))) * (tan(theta-pi/6.)/cos(theta
     1  -pi/6.)) * (one/(q*sin(three*theta)))
        dfthetads1 = dfthetads * (third+cos(three*theta)*ds1/(two*q)
     1              -three*((ds1**2)+(sig4**2)+(sig6**2))/(two*(q**2)))
        dfthetads2 = dfthetads * (third+cos(three*theta)*ds2/(two*q)
     1              -three*((sig4**2)+(ds2**2)+(sig5**2))/(two*(q**2)))
        dfthetads3 = dfthetads * (third+cos(three*theta)*ds3/(two*q)
     1              -three*((sig6**2)+(sig5**2)+(ds3**2))/(two*(q**2)))
        dfthetads4 = dfthetads * (      cos(three*theta)*sig4/(two*q)
     1         -three*((ds1*sig4)+(ds2*sig4)+(sig6*sig5))/(two*(q**2)))
        dfthetads5 = dfthetads * (      cos(three*theta)*sig5/(two*q)
     1         -three*((sig4*sig6)+(ds2*sig5)+(ds3*sig5))/(two*(q**2)))
        dfthetads6 = dfthetads * (      cos(three*theta)*sig6/(two*q)
     1         -three*((ds1*sig6)+(sig4*sig5)+(ds3*sig6))/(two*(q**2)))
C 
        fac1 = con1 * three/(two*q)
        fac2 = con1 * facT * facStrrate * facDctlDmgVol * sigmayield
C       Now finally calculating the tensor
        on1 = fac1*ds1
     1  - fac2 * (ftheta * dfitads1 + fita * dfthetads1)
C
        on2 = fac1*ds2
     1  - fac2 * (ftheta * dfitads2 + fita * dfthetads2)
C
        on3 = fac1*ds3
     1  - fac2 * (ftheta * dfitads3 + fita * dfthetads3)
C
        on4 = fac1*sig4
     1  - fac2 * (ftheta * dfitads4 + fita * dfthetads4)
C
        on5 = fac1*sig5
     1  - fac2 * (ftheta * dfitads5 + fita * dfthetads5)
C
        on6 = fac1*sig6
     1  - fac2 * (ftheta * dfitads6 + fita * dfthetads6)
C       Calculating the magnitude of the tensor
        onmag = sqrt(on1*on1+on2*on2+on3*on3+two*(on4*on4+on5*on5+on6
     1          *on6))
C
C 
C 
C NEWTON ITERATION
C       Preparing a variable for the equivalent plastic strain increment
        gammad  = zero
C       Preparing variables for the deviatoric stress tensor
        dds1    = ds1
        dds2    = ds2
        dds3    = ds3
        ssig4   = sig4
        ssig5   = sig5
        ssig6   = sig6
C       Preparing variables for the plastic strain increment tensor
        depsp1  = zero
        depsp2  = zero
        depsp3  = zero
        depsp4  = zero
        depsp5  = zero
        depsp6  = zero
C       Preparing a variable for the radius of elastic domain
        radiuss = radius * con1
C       Checking if damaging
        facdam = zero
        if (stateOld(i,4).gt.(one+half)
     1  .and. stateOld(i,4).lt.(two+half)) then
          facdam = one
        endif
        charLengthh = charLength(i)
        damold      = stateOld(i,2)
        sonset      = stateOld(i,2)
C       Newton iteration starts only if actually yielding
        IF (facyld.GT.half) THEN
        DO 130 kewton=1,newton
           rhs     = dsmag - radiuss
           dgammad = rhs/(zero-(twomu/dsmag)*(on1*dds1+on2*dds2+on3
     1               *dds3
     1     +two*(on4*ssig4+on5*ssig5+on6*ssig6))
     1     - (hard*onmag*fita*ftheta) + twomu*sigmayield*(
     1     ftheta*(on1*dfitads1+on2*dfitads2+on3*dfitads3
     1     +two*(on4*dfitads4+on5*dfitads5+on6*dfitads6)) +
     1     fita*(on1*dfthetads1+on2*dfthetads2+on3*dfthetads3
     1     +two*(on4*dfthetads4+on5*dfthetads5+on6*dfthetads6))
     1     ))
           gammad  = gammad - dgammad
           depsp1  = gammad * on1
           depsp2  = gammad * on2
           depsp3  = gammad * on3
           depsp4  = gammad * on4
           depsp5  = gammad * on5
           depsp6  = gammad * on6
C
C          DAMAGE
           dam = damold + facdam*(charLengthh*stateOld(i,3)
     1                     /(two*gf)) * gammad
           facDctlDmgDev = one - dam
           facDctlDmgVol = facDctlDmgDev
C
C          Updating flow stress
           call ahard(sigmayield,hard,stateOld(i,1)+gammad,props(33)
     1               ,nvalue)
C
C          Updating the stress tensor inside the newton iteration
           dds1    = ds1  - facDctlDmgDev * twomu * depsp1
           dds2    = ds2  - facDctlDmgDev * twomu * depsp2
           dds3    = ds3  - facDctlDmgDev * twomu * depsp3
           ssig4   = sig4 - facDctlDmgDev * twomu * depsp4
           ssig5   = sig5 - facDctlDmgDev * twomu * depsp5
           ssig6   = sig6 - facDctlDmgDev * twomu * depsp6
           dsmag = sqrt( dds1**2 + dds2**2 + dds3**2 + two*ssig4**2
     1     + two*ssig5**2
     1     + two*ssig6**2 )
C 
C Eq. (2) - Calculating the 2nd invariant of the stress tensor
           q   = dsmag / con1
C
C Eq. (3) - Calculating the 3rd invariant of the stress tensor
           r   = ( three/twoThirds * (dds1*(dds1**2 + ssig4**2 + ssig6
     1     **2) + 2.
     1     *ssig4*(dds1*ssig4 + dds2*ssig4 + ssig6*ssig5) + 2.*ssig6
     1     *(dds1*ssig6 + ssig4
     1     *ssig5 + dds3*ssig6) + dds2*(ssig4**2 + dds2**2 + ssig5**2)
     1     + 2.*ssig5
     1     *(ssig4*ssig6 + dds2*ssig5 + dds3*ssig5) + dds3*(ssig6**2
     1     + ssig5**2
     1     + dds3**2)) )**third
C
C Eq. (5) - Calculating the dimensionless hydrostatic pressure eta
           sigmean = third * ( dds1 + dds2 + dds3 ) + facDctlDmgVol
     1               * akappa * epsilontrace
           eta = sigmean / q
C
C Eq. (6) - Calculating the Lode angle theta
           cosine=(r/q)**3
C          Assuring that -1<cosine<1
           if (cosine.gt.one) cosine=one
           if (cosine.lt.(zero-one)) cosine=zero-one
           theta = third*acos(cosine)
C
C Eq. (7) - Calculating the normalized Lode angle thetabar
           thetabar = one - 6.*theta/pi
C
C Eq. (12) - Calculating the Lode angle dependant parameter gamma
           gamma = con2 * (one/cos(theta-pi/6.) - one)
C          Eq. (11) - Calculating radius of elastic zone
           fita   = one - ceta*(eta-eta0)
           ftheta = cthetas + ((cthetaax-cthetas)* (gamma- ((gamma**(om
     1     +one))/(om+one)) ) )
           radiuss = con1 * sigmayield * fita * ftheta * facT
     1     * facStrrate
     1     * facDctlDmgVol
C 
C Eq. (20) - Calculating the tensor of the normal direction of plastic flow
C          Avoiding a divide by zero by avoiding that theta = zero			 
           if (theta.eq.zero) theta = theta + tol
C          Derivative of triaxiality correction function w.r.t. stress
           dfitads  = ceta * three*eta/(two*(q**2))
           dfitads1 = dfitads * dds1
           dfitads2 = dfitads * dds2
           dfitads3 = dfitads * dds3
           dfitads4 = dfitads * ssig4
           dfitads5 = dfitads * ssig5
           dfitads6 = dfitads * ssig6
C          Derivative of Lode angle correction function w.r.t. stress
          dfthetads  = (cthetaax-cthetas) * (one -(gamma**om)) * (three
     1    *sqrt(three)/(two-sqrt(three))) * (tan(theta-pi/6.)/cos(theta
     1    -pi/6.)) * (one/(q*sin(three*theta)))
          dfthetads1 = dfthetads * (third+cos(three*theta)*dds1/(two*q)
     1           -three*((dds1**2)+(ssig4**2)+(ssig6**2))/(two*(q**2)))
          dfthetads2 = dfthetads * (third+cos(three*theta)*dds2/(two*q)
     1           -three*((ssig4**2)+(dds2**2)+(ssig5**2))/(two*(q**2)))
          dfthetads3 = dfthetads * (third+cos(three*theta)*dds3/(two*q)
     1           -three*((ssig6**2)+(ssig5**2)+(dds3**2))/(two*(q**2)))
         dfthetads4 = dfthetads * (      cos(three*theta)*ssig4/(two*q)
     1   -three*((dds1*ssig4)+(dds2*ssig4)+(ssig6*ssig5))/(two*(q**2)))
         dfthetads5 = dfthetads * (      cos(three*theta)*ssig5/(two*q)
     1   -three*((ssig4*ssig6)+(dds2*ssig5)+(dds3*ssig5))/(two*(q**2)))
         dfthetads6 = dfthetads * (      cos(three*theta)*ssig6/(two*q)
     1   -three*((dds1*ssig6)+(ssig4*ssig5)+(dds3*ssig6))/(two*(q**2)))
C         
           fac1 = con1 * three/(two*q)
           fac2 = con1 * facT * facStrrate * facDctlDmgVol * sigmayield
C          Now finally calculating the tensor
           on1 = fac1*dds1
     1     - fac2 * (ftheta * dfitads1 + fita * dfthetads1)
C         
           on2 = fac1*dds2
     1     - fac2 * (ftheta * dfitads2 + fita * dfthetads2)
C         
           on3 = fac1*dds3
     1     - fac2 * (ftheta * dfitads3 + fita * dfthetads3)
C         
           on4 = fac1*ssig4
     1     - fac2 * (ftheta * dfitads4 + fita * dfthetads4)
C         
           on5 = fac1*ssig5
     1     - fac2 * (ftheta * dfitads5 + fita * dfthetads5)
C         
           on6 = fac1*ssig6
     1     - fac2 * (ftheta * dfitads6 + fita * dfthetads6)
C          Calculating the magnitude of the tensor
           onmag = sqrt(on1*on1+on2*on2+on3*on3+two*(on4*on4+on5*on5
     1             +on6*on6))
C 
           if(abs(rhs).lt.tol*radius) goto 140
 130    CONTINUE
        WRITE(6,2) newton
 2      FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     1      'CONVERGE AFTER ',I3,' ITERATIONS')
 140    CONTINUE
        ENDIF
C 
        radiuss = radiuss/con1
        stateNew(i,2) = dam
C 
C 
C 
C 
C Storing equivalent plastic strain increment
        stateNew(i,17) = gammad
C
C Updating and storing equivalent plastic strain
        stateNew(i,1) = stateOld(i,1) + gammad
C
C Updating and storing the plastic strain tensor (states 11-16)
        stateNew(i,11) = stateOld(i,11) + depsp1 /twoThirds
        stateNew(i,12) = stateOld(i,12) + depsp2 /twoThirds
        stateNew(i,13) = stateOld(i,13) + depsp3 /twoThirds
        stateNew(i,14) = stateOld(i,14) + depsp4 /twoThirds
        stateNew(i,15) = stateOld(i,15) + depsp5 /twoThirds
        stateNew(i,16) = stateOld(i,16) + depsp6 /twoThirds
C
C Updating and storing the equivalent plastic strain rate
        stateNew(i,22) = stateOld(i,22) + one
        stateNew(i,23) = stateOld(i,23) + gammad
        stateNew(i,24) = stateOld(i,24) + dt
        if (strrInc.gt.half .and. strrInc.lt.(one+half)) then
          stateNew(i,21) = gammad/dt
          stateNew(i,27) = gammad
        else
          stateNew(i,21) = (stateNew(i,23)/stateNew(i,24)
     1                   +  stateOld(i,21))/two
          stateNew(i,27) = (stateNew(i,23)/stateNew(i,22)
     1                   +  stateOld(i,27))/two
          if (stateNew(i,22).gt.strrInc-half .and. stateNew(i,22)
     1    .lt.strrInc+half) then
            stateNew(i,22) = zero
            stateNew(i,23) = zero
            stateNew(i,24) = zero
          endif
        endif
        stateNew(i,20) = gammad / dt
C
C Updating temperature
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
          tempNew(i) = t0
     1    +(eta2/(density(i)*cp))*radiuss*gammad/facDctlDmgDev
        else
          tempNew(i) = tempOld(i)
     1    +(eta2/(density(i)*cp))*radiuss*gammad/facDctlDmgDev
        endif
C
C Calculating equivalent plastic strain to failure
        epsilonf = (half * ( d1*exp(zero-d2*eta) + d5*exp(zero-d6*eta))
     1                     - d3*exp(zero-d4*eta) )
     1                     * thetabar**2
     1           +  half * ( d1*exp(zero-d2*eta) - d5*exp(zero-d6*eta))
     1                     * thetabar
     1           +  d3*exp(zero-d4*eta)
        stateNew(i,29) = epsilonf
C
C Ductile damage
        if (stateNew(i,1).ge.epsilonf .and. stateNew(i,1).gt.tol) then
          if (stateOld(i,3).le.tol .and. stateOld(i,3).ge.-tol) then
C           Registering the yield stress at the onset of damage for the first time
            stateNew(i,3) = radiuss
          else
C           Keeping the yield stress at the onset of damage for next steps
            stateNew(i,3) = stateOld(i,3)
          endif
C         Updating damage fraction variable
C          CHANGE
C          stateNew(i,2) = stateOld(i,2) + (charLength(i)*stateNew(i,3)
C     1                    /(two*gf)) * gammad
          stateNew(i,4) = two
        else
C
C         In case no damage is reached, yield stress at damage onset
C         and damage variable are kept for next step
          stateNew(i,3) = stateOld(i,3)
          stateNew(i,2) = stateOld(i,2)
        endif
C
C       Previous ductile damage
        if (stateOld(i,4).gt.(one+half)
     1  .and. stateOld(i,4).lt.(two+half)) then
          stateNew(i,4) = two
        endif
C
C Update stress
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * gammad / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
        sig1 = sig1 + fac1 - fac2 * on1 /con1
        sig2 = sig2 + fac1 - fac2 * on2 /con1
        sig3 = sig3 + fac1 - fac2 * on3 /con1
        sig4 = sig4        - fac2 * on4 /con1
        sig5 = sig5        - fac2 * on5 /con1
        sig6 = sig6        - fac2 * on6 /con1
C       Calculating invariants of stress tensor
        SI1 =   sig1 + sig2 + sig3
        SI2 =   sig1*sig2-sig4*sig4
     1        + sig1*sig3-sig6*sig6
     1        + sig2*sig3-sig5*sig5
        SI3 =   sig1*(sig2*sig3-sig5*sig5)
     1        - sig4*(sig4*sig3-sig5*sig6)
     1        + sig6*(sig4*sig5-sig2*sig6)
C       Preparing subvalues for calculating the principal stresses values
        cosine2 = (two*SI1*SI1*SI1-three*three*SI1*SI2+three*three
     1            *three*SI3)
     1            /(two*(SI1*SI1-three*SI2)**(threehalves))
C       Assuring that -1<cosine2<1
        if (cosine2.gt.one) cosine2=one
        if (cosine2.lt.(zero-one)) cosine2=zero-one
        alpha2 = acos(cosine2)
C       Calculating the principal stress values
        SP1 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three)
        SP2 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three+twoThirds*pi)
        SP3 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(twoThirds*pi-alpha2/three)
C       Fetching the highest of the principal stress values
        sigmamax = max(abs(SP1),abs(SP2),abs(SP3))
        stateNew(i,26) = sigmamax
C       Calculating new softening correction terms due to ductile damage
        stateNew(i,28)= one
        if (stateNew(i,2).ge.Dcrit) then
C          DAMAGE ON
          facDctlDmgDev = zero
          facDctlDmgVol = zero
C          facDctlDmgDev = one
C          facDctlDmgVol = one
          stateNew(i,2) = Dcrit
C          DAMAGE ON
          stateNew(i,28)= zero
        else
C          CHANGE
C          facDctlDmgDev = one - stateNew(i,2)
C          facDctlDmgVol = facDctlDmgDev
C          facDctlDmgDev = one
C          facDctlDmgVol = one
        endif
C       Elements under hydrostatic compression don't experience spherical damage
        if (sig1+sig2+sig3.lt.zero .and. stateNew(i,2).lt.Dcrit) then
          facDctlDmgVol = one
        endif
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * gammad / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
        stressNew(i,1) = sig1 + fac1 - fac2 * on1 /con1
        stressNew(i,2) = sig2 + fac1 - fac2 * on2 /con1
        stressNew(i,3) = sig3 + fac1 - fac2 * on3 /con1
        stressNew(i,4) = sig4        - fac2 * on4 /con1
        stressNew(i,5) = sig5        - fac2 * on5 /con1
        stressNew(i,6) = sig6        - fac2 * on6 /con1
C       Brittle damage
C       DAMAGE ON
        if (sigmamax.ge.sigdmg .and. stressNew(i,1).gt.zero) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4)  = one
          stateNew(i,28) = zero
        endif
C
C
  100 continue
C
      return
      end
C
C
      subroutine ahard(sigmayield,hard,eqplas,table,nvalue)
C
      include 'vaba_param.inc'
      dimension table(2,nvalue)
C
C     Set yield stress to last value of table, hardening to zero
      sigmayield=table(1,nvalue)
      hard=0.0
C
C     If more than one entry, search table
      if(nvalue.gt.1) then
        do 10 k1=1,nvalue-1
          eqpl1=table(2,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(2,k1)
            if(eqpl1.le.eqpl0) then
              write(6,7)
 7            format(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     1               'ENTERED IN ASCENDING ORDER,')
C
C             Subroutine XIT terminates execution and closes all files
              call XIT
            endif
            deqpl=eqpl1-eqpl0
            sigmayield0=table(1,k1)
            sigmayield1=table(1,k1+1)
            dsigmayield=sigmayield1-sigmayield0
            hard=dsigmayield/deqpl
            sigmayield=sigmayield0+(eqplas-eqpl0)*hard
            goto 20
          endif
 10     continue
 20     continue
        if(eqplas.gt.table(2,nvalue)) then
          hard=(table(1,nvalue)-table(1,nvalue-1))
     1        /(table(2,nvalue)-table(2,nvalue-1))
          sigmayield=table(1,nvalue)+(eqplas-table(2,nvalue))*hard
        endif
      endif
      return
C
C Iteration ends here
      end
