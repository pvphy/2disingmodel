 	program isingmodel
	implicit none
	double precision::T,dt,E,random1,ptrial, Estart,Eflip,Efinal,de,E_T,Et2,m,avg_E,avg2_E,Cv,mm,ef,e1,e2,avg_mag
	!i=col,j=row
	integer :: i,j,rt,lt,up,dw
	integer :: nmc,nn,it,nnn
	integer, parameter :: J1 = -1  
	integer,dimension(:,:),allocatable :: spin(:,:) 
	print*,'give sites'
	read*,nn
	print*,'give mc step'
	read*,nnn
	allocate(spin(nn,nn))
	spin=0
	
	do i=1,nn
	 do j=1,nn
      call random_number(random1)
      if (random1 .lt. 0.5) then
       spin(i,j) = +1
       else 
       spin(i,j) = -1
      end if
	 end do
	end do
	!write(*,*)spin(1,1)
	! T=10.010d0
	
	!do T=5,0.1,-0.05
    T=2.010d0
	
	do iT=1,30
	t=t+dt
	
	if(iT.le.10) T=T-0.1
	if ((iT.gt.10).and.(iT.le.20)) T=T-0.09
	if ((iT.gt.20).and.(iT.le.30)) T=T-0.009
	
	e1=0.0d0
	e2=0.0d0
	m=0.d0
	!mm=0.d0
	
	 do nmc=1,nnn
	 e_t=0.0d0
      do i=1,nn
       do j=1,nn
            rt = i+1
            lt = i-1
            dw= j-1
            up= j+1

            if (i.eq.1) lt= nn
            if (i.eq.nn) rt = 1
            if (j.eq.1) dw= nn
            if (j.eq.nn) up= 1
	     	Estart = J1*spin(i,j)*(spin(rt,j)+spin(lt,j)+spin(i,dw)+spin(i,up))
	        Eflip=J1*(-spin(i,j))*(spin(rt,j)+spin(lt,j)+spin(i,dw)+spin(i,up))
            de=Eflip-Estart
        

	      if(de.le.0) then 
           spin(i,j) = -spin(i,j)
           Efinal=Eflip
           else  
           call random_number(random1)
           ptrial = EXP(-(de)/T)

          if (random1 .lt. ptrial) then
           spin(i,j) = -spin(i,j)
           Efinal=Eflip
           else
           Efinal=Estart
          end if
           E_T=E_T+Efinal 
           Et2=Et2+ Efinal**2
         
	    end if 	    
        end do
        enddo
        if(nmc.gt.(nnn/2))then
        do i=1,nn
        do j=1,nn   
        m=m+spin(i,j) 
        enddo        
        enddo
     do i=1,nn
       do j=1,nn
            rt = i+1
            lt = i-1
            dw= j-1
            up= j+1

            if (i.eq.1) lt= nn
            if (i.eq.nn) rt = 1
            if (j.eq.1) dw= nn
            if (j.eq.nn) up= 1
	     	Ef = J1*spin(i,j)*(spin(rt,j)+spin(lt,j)+spin(i,dw)+spin(i,up))
           e1=e1+(ef/2.0)
	       e2=e2+((ef/2.0)**2)
           enddo
           enddo
      endif
      enddo	
	!write(30,*) t,abs(m)/dble(n*n*(nn/2)),avg_e/dble(n*n*(nn/2))
            avg_mag=abs(m)/dble(nn*nn*(nnn/2))
            avg_E=e1/dble(nn*nn*(nnn/2)) 
            avg2_E=e2/dble(nn*nn*(nnn/2)) 
            cv=(avg2_e-(avg_e)**2)/((t**2))
            write(11,*) t,avg_mag,avg_e,cv	      
	end do   
	
	end program isingmodel
