!---------------------------------
! define the initial location
!--------------------------------
  implicit none
  double precision,intent(in)::R_ta,R_im
  double precision,intent(inout):: xinit1,yinit1,vxinit1,vyinit1
  double precision,intent(inout):: xinit2,yinit2,vxinit2,vyinit2
  double precision :: m_ta,m_im,ratio,guzai,m_tot!L_in must be in insph
  double precision :: L_tochu,GGG
  double precision :: xxx(3,3),vvv(3,3),dtt,ttt,aa(4),rr(2),dis
  double precision :: rg(2),m1,m2,myu,fr,vvr(2,2),min
  double precision :: maxvel,revx,revy,rabs,vabs,vesc,velmugen,velhere
  double precision :: RE,v_im,v_ta,arg,ME,Rearth
  m_ta=mas(1)*npart_l
  m_im=mas(npart_l+1)*(npart-npart_l)
  m_tot=m_ta+m_im
  ratio=m_im/m_tot
  ME=6.d0*10.d0**24.d0
  Rearth=6.4d0*10.d0**6.d0
  GGG=6.670d0*10.d0**(-11.d0)
  v_im=vel_inpt*dsqrt(2.d0*GGG*m_tot/(R_ta+R_im))*(1.d0-ratio)
  v_ta=vel_inpt*dsqrt(2.d0*GGG*m_tot/(R_ta+R_im))*ratio
  arg=dasin(scaled_b)
 print *,v_ta/vel_inpt/ratio
  write(*,*) 'total mass=',m_tot/ME,'target_rad=',R_ta/Rearth
  write(*,*) 'impactor_rad=',R_im/Rearth,'mass ratio=',ratio
  write(*,*) 'v_imp/v_esc=',vel_inpt,'b',scaled_b
  write(*,*)'vesc=',dsqrt(2.d0*GGG*m_tot/(R_ta+R_im))
  write(41,*)'vesc=',dsqrt(2.d0*GGG*m_tot/(R_ta+R_im))
  close(41)
  !stop
  m1=m_ta
  m2=m_im
  maxvel=0.d0
  xxx(1,1)=0.d0!-R_ta  ! target x position
  xxx(1,2)=0.d0  ! target y position
  xxx(2,1)=(R_ta+R_im)*dcos(arg)  ! impactor x position
  xxx(2,2)=(R_ta+R_im)*dsin(arg)  ! impactor y position
     dis=(xxx(1,1)-xxx(2,1))**2.d0
     dis=dis+(xxx(1,2)-xxx(2,2))**2.d0
     dis=dsqrt(dis)
  vvv(1,1)=v_ta!*cos(arg)  ! target x velocity
  vvv(1,2)=0.d0!v_ta*sin(arg)  ! target y velocity
  vvv(2,1)=-v_im!*cos(arg)  ! impactor x velcocity
  vvv(2,2)=0.d0!-v_im*sin(arg)  ! impactor y velocity
  ttt=0.d0  ! time
  dtt=-0.1d0  !timestep
  RE=R_ta  ! radius earth = radius target
  do while(ttt>-5000.d0)!dis>1.d0*RE)
     dis=(xxx(1,1)-xxx(2,1))**2.d0
     dis=dis+(xxx(1,2)-xxx(2,2))**2.d0
     dis=dsqrt(dis)  ! calculate the distance between target and impactor
     revx=xxx(1,1)-xxx(2,1)  ! calculate the x distance between target and impactor
     revy=xxx(1,2)-xxx(2,2)  ! calculate the y distance between target and impactor
     if(dis<=min)then
        min=dis
     endif
     if(dis/RE>=3.0d0)then
        xinit1=xxx(1,1)/RE
        yinit1=xxx(1,2)/RE
    xinit2=xxx(2,1)/RE
    yinit2=xxx(2,2)/RE
        vxinit1=vvv(1,1)
        vyinit1=vvv(1,2)
        vxinit2=vvv(2,1)
        vyinit2=vvv(2,2)
    print *,xinit1,yinit1,xinit2,yinit2,"place"
        return
     endif
     aa(1)=dtt*(GGG*m2/dis/dis/dis*(xxx(2,1)-xxx(1,1)))
     aa(2)=dtt*(GGG*m2/dis/dis/dis*(xxx(2,2)-xxx(1,2)))
     aa(3)=-aa(1)*m1/m2
     aa(4)=-aa(2)*m1/m2
     vvv(1,1)=vvv(1,1)+aa(1)
     vvv(1,2)=vvv(1,2)+aa(2)
     vvv(2,1)=vvv(2,1)+aa(3)
     vvv(2,2)=vvv(2,2)+aa(4)
     xxx(1,1)=xxx(1,1)+dtt*vvv(1,1)
     xxx(1,2)=xxx(1,2)+dtt*vvv(1,2)
     xxx(2,1)=xxx(2,1)+dtt*vvv(2,1)
     xxx(2,2)=xxx(2,2)+dtt*vvv(2,2)
     ttt=ttt+dtt
  enddo
  end subroutine define_place