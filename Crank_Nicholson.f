C       MARTA XIULAN ARIBÓ HERRERA
C       DESCRIPCIÓ: Resoldrem l'equació de calor
c-----------------------------------------------------------------------
        PROGRAM P10
        IMPLICIT NONE
c--------------------------DECLARACIÓ DE VARIABLES----------------------
        integer Nx, k, kmax, i,C,Nt,j
        parameter (Nx = 50, kmax = 1000, Nt=6000)
        double precision L,h,tx1(0:Nx),tx0(0:Nx),t_a,beta(3),alpha_f,dx
        double precision tx2(0:Nx),w,dt,temps,alpha_v(2),alpha_au
        parameter (L = 150.d0, h = L/Nx, t_a = 22.d0,dx = L/Nx )!cm
        parameter (temps = 4500, dt =temps/Nt)
        parameter (alpha_f =2.2d0*10**(-1.d0))
        double precision error(0:Nx)
        double precision u0(0:Nx),u1(0:Nx),Bu(1:Nx-1)
        double precision BB(1:Nx-1,1:Nx-1),AP(1:Nx-1),A0(1:Nx-1),
     *AM(1:Nx-1),r,b, aux(1:Nx-1)
        beta = (/0.00004d0,0.0003d0, 0.0025d0/) !Diferents valors de beta
        open(1,file ="check1.dat") !fitxers auxiliars
        open(2,file ="check2.dat")
        open(3,file ="check3.dat")
        open(4,file ="check4.dat")
        !condicions de contorn
        tx0(Nx)= 280.d0-t_a !ºC !condicions inicials
        tx0(0)= 0.d0
        print*,dt,"ashe",h
c        print*, tx0
c----------------------------PRIMERA PART-------------------------------
C       Calculem el perfil de temperatures per diferents valors de beta
c       en estat estacionari
c-----------------------------------------------------------------------
        DO C= 1,3
            write(1,*)"#Beta",beta(c)
            do i = 1,Nx-1
            tx0(i) = 150.d0
            enddo
            do k = 1, kmax
            tx1 = tx0
            do i = 1, Nx-1!!!
                w = ((beta(c)/alpha_f)*(h**2.d0))+2.d0
                tx0(i)=(tx0(i+1)+tx0(i-1))*(1.d0/w)
            enddo!!!
            error = abs(tx1-tx0)
            if (maxval(error).LE.0.0000001d0) then
                print*, "s'ha arribat a la convergència desitjada",
     *maxval(error)
                exit
            endif
            enddo 

            do i = 0,Nx
                 write(1,*) i,tx0(i)+t_a
            enddo 
            write(1,"(a1)")
            write(1,"(a1)")
        ENDDO 
        call system ("gnuplot -p plot1.gnu")
c---------------------------- SEGONA PART A -------------------------------
C       CRANCK NICHOLSON Au1 =Bu0
c-----------------------------------------------------------------------
        !construim les matrius A i B
        write(2,*) "#APARTAT 2A"
        r = alpha_f*dt/(h**2.d0)
        b = 0.00004d0
        BB = 0.D0!tot a zero menys la diagonal
        !Construcció de les matrius tridiagonals, i les tres diagonals en forma
        !de vectors per poder introduir a la subrutina tridiag.
        do i=1,Nx-1 
            do j = 1,Nx-1
            	if (i.EQ.j) then
                BB(i,j)= 1.d0-(b*dt/2.d0)-r
                endif
                if (i.EQ.(j+1)) then 
                    BB(i,j)= +r/2.d0
                endif
                if ((i+1).EQ.j) then 
                    BB(i,j)= +r/2.d0
                endif
            !Matrius diagonals per posarles a la subrutina tridiagonal
                    A0(i)= 1.d0+(b*dt/2.d0)+r
                    AP(i)= -r/2.d0
                    AM(i)= -r/2.d0
             enddo 
        enddo 
        u0 = 0.d0 
        u0(0) = 0.d0 !condicions de contorn fixades
        u0(Nx) = 280.d0-T_a !condicions de contorn fixades
        u1 = u0

        AP(Nx-1) = 0.d0 !extrems de les matrius
        AM(1) = 0.d0 

        DO K = 1, Nt
            do i = 1, Nx-1
            Bu(i) = 0.d0
            do j = 1,Nx-1
               Bu(i)= Bu(i)+(BB(i,j)*u0(j)) !multiplicació de A*U1
            enddo
            enddo 
            Bu(1) =Bu(1)+r*u0(0)
            Bu(Nx-1)=Bu(Nx-1)+r*u0(Nx)
            !calculamos la operación A·U1 = B·U0 con la subrutina tridiag
            do i = 1,Nx-1 !vector auxiliar per passar a la tridiag els vectors entremig de u1
                aux(i)= u1(i) !i que les dimensions quadrin
            enddo

            call tridiag(AM,A0,AP,Bu,aux,Nx-1)

            do i = 1,Nx-1
                u1(i)= aux(i) !es retornen els valors a u1
            enddo
            u0 = u1 !es machaquen els vectors

            write(2,*)k*dt,i*h,u1(2)+T_A,u1(14)+T_A,u1(42)+t_a !escriptura de dades
            !en el fitxer

        ENDDO 

        call system ("gnuplot -p plot2.gnu") !evolució temporal pels punts

c---------------------------- SEGONA PART B -------------------------------
C       Es calcula l'evolució temporal mitjana per cada temps per diferents
c       alpha's.
c-----------------------------------------------------------------------
        write(3,*) "#APARTAT 2B"
        alpha_au=1.29d0 !
        alpha_v=(/alpha_f,alpha_au/) !per compactar codi es posen en un vector
        do c = 1,2
            r = alpha_v(c)*dt/(h**2.d0)
            b = 0.00015d0
            BB = 0.D0!tot a zero menys la diagonal
            print*,2.d0*(1.d0+(b*Nt/2.d0)+r), r
            do i=1,Nx-1
                do j = 1,Nx-1
                    if (i.EQ.j) then
                    BB(i,j)= 1.d0-(b*dt/2.d0)-r
                    endif
                    if (i.EQ.(j+1)) then 
                        BB(i,j)= +r/2.d0
                    endif
                    if ((i+1).EQ.j) then 
                        BB(i,j)= +r/2.d0
                    endif
                !Matrius diagonals per posarles a la subrutina tridiagonal
                        A0(i)= 1.d0+(b*dt/2.d0)+r
                        AP(i)= -r/2.d0
                        AM(i)= -r/2.d0
                 enddo 
            enddo 
            u0 = 0.d0
            u0(0) = 0.d0
            u0(Nx) = 280.d0-T_a
            u1 = u0

            AP(Nx-1) = 0.d0
            AM(1) = 0.d0 
            DO K = 1, Nt
                do i = 1, Nx-1
                Bu(i) = 0.d0
                do j = 1,Nx-1
                   Bu(i)= Bu(i)+(BB(i,j)*u0(j))
                enddo
                enddo 
                Bu(1) =Bu(1)+r*u0(0)
                Bu(Nx-1)=Bu(Nx-1)+r*u0(Nx)
                !calculamos la operación A·U1 = B·U0 con la subrutina tridiag
                do i = 1,Nx-1
                    aux(i)= u1(i)
                enddo

                call tridiag(AM,A0,AP,Bu,aux,Nx-1)

                do i = 1,Nx-1
                    u1(i)= aux(i)
                enddo
                u0 = u1

                write(3,*)k*dt, sum(u1+T_a)/L,c
            
            ENDDO 
            write(3,"(a1)")
            write(3,"(a1)")
        enddo
        call system ("gnuplot -p plot3.gnu") !gràfica de l'evolució temporal de temperatures mitjana
c        print*,"je",u1
c---------------------------- SEGONA PART C -------------------------------
C       Evolució temporal de la temperatura per una beta donada 0.002
c-----------------------------------------------------------------------
        !construim les matrius A i B
        write(2,*) "#APARTAT 2A"
        r = alpha_au*dt/(h**2.d0)
        b = 0.002d0
        BB = 0.D0!tot a zero menys la diagonal
c        print*,2.d0*(1.d0+(b*Nt/2.d0)+r), r
        do i=1,Nx-1
            do j = 1,Nx-1
                if (i.EQ.j) then
                BB(i,j)= 1.d0-(b*dt/2.d0)-r
                endif
                if (i.EQ.(j+1)) then 
                    BB(i,j)= +r/2.d0
                endif
                if ((i+1).EQ.j) then 
                    BB(i,j)= +r/2.d0
                endif
            !Matrius diagonals per posarles a la subrutina tridiagonal
                    A0(i)= 1.d0+(b*dt/2.d0)+r
                    AP(i)= -r/2.d0
                    AM(i)= -r/2.d0
             enddo 
        enddo 
        u0 = 0.d0
        u0(0) = 0.d0
        u0(Nx) = 280.d0-T_a
        u1 = u0

        AP(Nx-1) = 0.d0
        AM(1) = 0.d0 

        DO K = 1, Nt
            do i = 1, Nx-1
            Bu(i) = 0.d0
            do j = 1,Nx-1
               Bu(i)= Bu(i)+(BB(i,j)*u0(j))
            enddo
            enddo 
            Bu(1) =Bu(1)+r*u0(0)
            Bu(Nx-1)=Bu(Nx-1)+r*u0(Nx)
            !calculamos la operación A·U1 = B·U0 con la subrutina tridiag
            do i = 1,Nx-1
                aux(i)= u1(i)
            enddo

            call tridiag(AM,A0,AP,Bu,aux,Nx-1)

            do i = 1,Nx-1
                u1(i)= aux(i)
            enddo
            u0 = u1
            do i = 0,Nx
            write(4,*)dt*k,i*h,u1(i)+T_a
            enddo             
            write(4,"(a1)")
            write(4,"(a1)")
        ENDDO 
        call system("gnuplot -p plot4.gnu")
C-----------------------------------------------------------------------
        END PROGRAM

c-------------------------Subrutines i funcions-------------------------

C--------------------- SUBRUTINA TRIDIAGONAL -----------------------
c 
c Solves the problem T psi =R
c 
c where T is a tridiagonal matrix, A (lower), B (central), C (upper)
c  A   0 a1, a2, a3, ...., aIMAX
c  B   b1 b2, b3, b4, ...., bIMAX
c  C   c1, c2, c3, c4, ....,0 

        SUBROUTINE TRIDIAG(A,B,C,R,PSI,IMAX)
        IMPLICIT double precision (A-H,K,O-Z)
        IMPLICIT INTEGER (I-J , L-N)
        double precision  BET
        double precision  GAM(4001)
        double precision A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

        IF(B(1).EQ.0.) PAUSE
        BET=B(1)
        PSI(1)=R(1)/BET
        DO 11 J=2,IMAX
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0) PAUSE
        PSI(J)=(R(J)-A(J)*PSI(J-1))/BET
11      CONTINUE

        DO 12 J=IMAX-1,1,-1
        PSI(J)=PSI(J)-GAM(J+1)*PSI(J+1)
12      CONTINUE

       RETURN
       END




