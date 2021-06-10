do i = 1, npart
    if(mas(i)<1.d0)then !a part of a group
        continue
    else
        rad(i) = dsqrt(rx(i) * rx(i) + ry(i) * ry(i) + rz(i) * rz(i))! radial distance from the center
        vel(i) = dsqrt(vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i)) ! velocity
        hhh(i) = dsqrt((ry(i) * vz(i) - rz(i) * vy(i))**2.d0 + (rz(i) * vx(i) - rx(i) * vz(i))**2.d0 + (rx(i) * vy(i) - ry(i) * vx(i))**2.d0) ! angular momentum
        Lz2(i) = (rx(i) * vy(i) - ry(i) * vx(i))! angular momentum in the z-direction
        inc(i) = dacos(dabs(Lz2(i) / hhh(i))) !inclination
        rad_av(i) = Lz2(i) * Lz2(i) / GG / Mpe
        if((rad(i)<=aa).or.(rad_av(i))<=aa)then
            Mpe2 = Mpe2 + mas(i)! the particle becomes part of the planet
            Lpla2 = Lpla2 + mas(i) * (rx(i) * vy(i) - ry(i) * vx(i)) ! angular momentum of the planet
            vrr = (-vx(i) * ry(i) + rx(i) * vy(i)) / dsqrt(rx(i) * rx(i) + ry(i) * ry(i))
        else if(type_here(i)==0)then !the particle belongs to the planet. Here, type_here(i)==- means that it was previously indentifed that the particle is part of the planet
            Mpe2 = Mpe2 + mas(i)
            Lpla2 = Lpla2 + mas(i) * (rx(i) * vy(i) - ry(i) * vx(i))
            vrr = (-vx(i) * ry(i) + rx(i) * vy(i)) / dsqrt(rx(i) * rx(i) + ry(i) * ry(i))
        else
            ene(i) = 0.5d0 * mas(i) * vel(i) * vel(i) - GG * Mpe * mas(i) / rad(i) ! calculating the orbital energy
            !eee(i)=dsqrt(1.d0+2.d0*ene(i)*hhh(i)*hhh(i)/mas(i)/GG/GG/Mpe/Mpe)
            eee(i) = dsqrt(1.d0 + 2.d0 * ene(i) * Lz2(i)**2.d0 / mas(i) / GG / GG / Mpe / Mpe) ! calculating the eccentricity
            aaa(i) = -GG * Mpe * mas(i) / 2.d0 / ene(i) ! calculating semi-major axis
            rad_org(i) = aaa(i)
            if(eee(i)<1.d0)then ! eliptic
                qqq(i) = aaa(i) * (1.d0 - eee(i))
                if(qqq(i)<=aa)then !a particle will fall onto the planet at some point
                    Mpe2 = Mpe2 + mas(i)
                    Lpla2 = Lpla2 + mas(i) * (rx(i) * vy(i) - ry(i) * vx(i))
                    vrr = (-vx(i) * ry(i) + rx(i) * vy(i)) / dsqrt(rx(i) * rx(i) + ry(i) * ry(i))
                else !disk particle
                    aaa(i) = rad_org(i) * (1.d0 - eee(i)**2.d0) * dcos(inc(i))**2.d0 ! this is assuming that eccentricity and angular momentum are damped
                    internal = internal + mas(i) * uu(i)
                    Mdisk = Mdisk + mas(i)
                    Ldisk = Ldisk + mas(i) * (rx(i) * vy(i) - ry(i) * vx(i))
                    disk(m) = i
                    m = m + 1
                endif
            else  ! parabolic orbit
                Mesc = Mesc + mas(i)
                Lesc = Lesc + mas(i) * (rx(i) * vy(i) - ry(i) * vx(i))
                escpat(j) = i
                j = j + 1
            endif
        endif
    endif
enddo
II = 2.d0 / 5.d0 * Mpe2 * aa * aa !morment of inertia
omega = Lpla2 / II !angular velocity of the planet
kepler = dsqrt(GG * Mpe2 / aa / aa / aa) !Keplerian angular velocity at r=aa
ff2 = 2.5d0 * (omega / kepler)**2.d0 / (1.d0 + (2.5d0 - 1.5d0)**2.d0)
aa_old = aa
aa = (3.d0 / 4.d0 / 3.141592d0 * Mpe2 / rho_ct / (1.d0 - ff2))**(0.33333d0) ! new aa
Mpe = Mpe2
enddo