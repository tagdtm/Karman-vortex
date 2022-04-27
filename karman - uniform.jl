using Printf

# only for ununiformed mesh
function mesh_y(ny,rx,ALG)   
    
    AT = Float64(log((1+ALG)/(1-ALG))/2)
    y1 = Int64((ny-rx)/2)
    y2 = Int64(117)
    y3 = Int64((ny-y2+1)/2)
    yp = zeros(Float64,ny)
    dy = zeros(Float64,ny)
    dya = zeros(Float64,ny)

    for i = 1:y3
        ETA = AT*(-1+2*(i+y1)/(2*y1))
        yp[i] = tanh(ETA)/ALG*y1/ny*12.8
    end
        
    for i = y2+y3:ny
        ETA = AT*(-1+2*(i-y1-rx)/(2*y1))
        yp[i] = ((1+tanh(ETA)/ALG)*y1/ny+(y1+rx)/ny)*12.8
    end
        
    for i = 1+y3:y2+y3-1
        yp[i] = (yp[y3]+(i-y3)*(yp[y2+y3]-yp[y3])/y2)
    end
    
    for i = 1:ny
        if i ==1
            dy[i] = yp[i]
        else
            dy[i] = yp[i]-yp[i-1]
            dya[i] = (dy[i]+dy[i-1])/2
        end
    end
     
    return yp,dy,dya
end

function periodic(ny)
    jp = zeros(Int64,ny)
    jm = zeros(Int64,ny)
    for j = 1:ny
        if j == ny
            jp[j] = 1
            jm[j] = j-1
        elseif j == 1
            jp[j] = j+1
            jm[j] = ny
        else
            jp[j] = j+1
            jm[j] = j-1
        end 
    end
    return jp,jm
end

function convection(nx,ny,rx,dx,dy,dya,u,v,Con_x_p, Con_y_p,jp,jm) 
    Con_x_f = zeros(Float64,nx,ny)
    Con_y_f = zeros(Float64,nx,ny)
    JU1 = zeros(Float64,nx,ny)
    JU2 = zeros(Float64,nx,ny)
    JV1 = zeros(Float64,nx,ny)
    JV2 = zeros(Float64,nx,ny)
    u1 = zeros(Float64,nx,ny)
    u2 = zeros(Float64,nx,ny)
    v1 = zeros(Float64,nx,ny)
    v2 = zeros(Float64,nx,ny)
    A_x = zeros(Float64,nx,ny)
    A_y = zeros(Float64,nx,ny)
    
    for j = 1:ny
        for i = 1:nx
            Con_x_f[i,j] = Con_x_p[i,j]
            Con_y_f[i,j] = Con_y_p[i,j]
        end
    end
    
    for j = 1:ny
        for i = 1:nx
            if(i > 1)
                JU1[i,j] = dy[j]*(u[i-1,j]+u[i,j])/2
                u1[i,j] = (u[i-1,j] + u[i,j])/2
            end
            if(i < nx)
                JV1[i,j] = (dx[i]*v[i,j]+dx[i+1]*v[i+1,j])/2
                v2[i,j] = (v[i,j] + v[i+1,j])/2
            end
            JU2[i,j] = (dy[j]*u[i,j]+dy[jp[j]]*u[i,jp[j]])/2
            u2[i,j] = (u[i,j] + u[i,jp[j]])/2
            JV2[i,j] = dx[i]*(v[i,jm[j]]+v[i,j])/2
            v1[i,j] = (v[i,jm[j]] + v[i,j])/2
        end
    end
    
    for j = 1:ny
        for i = 2:nx-1
            if (x1-2 < i && i < x1+rx && y1-1 < j && j < y1+rx)
                # do nothing
            else
                J1 = dx[i]*dy[j]
                Con_x_p[i,j] = (-JU1[i,j]*u1[i,j]+JU1[i+1,j]*u1[i+1,j])/J1+(-JV1[i,jm[j]]*u2[i,jm[j]]+JV1[i,j]*u2[i,j])/J1
                A_x[i,j] = - (3*Con_x_p[i,j]-Con_x_f[i,j])/2
            end
            
            if (x1-1 < i && i < x1+rx && y1-2 < j && j < y1+rx)
                # do nothing
            else
                J2 = dx[i]*dya[jp[j]]
                Con_y_p[i,j] = (-JU2[i-1,j]*v1[i-1,j]+JU2[i,j]*v1[i,j])/J2+(-JV2[i,j]*v2[i,j]+JV2[i,jp[j]]*v2[i,jp[j]])/J2
                A_y[i,j] = - (3*Con_y_p[i,j]-Con_y_f[i,j])/2
            end
        end
    end
    
    return A_x, A_y, Con_x_p, Con_y_p 
end

function diffusion(nx,ny,rx,dx,dy,dya,Re,dt,u,v,A_x,A_y,p,x1,y1,jp,jm)
    qu = zeros(Float64,nx,ny)
    qv = zeros(Float64,nx,ny)
    sumqu = Float64(0)
    sumqv = Float64(0)
    Dif_x = Float64(0)
    Dif_y = Float64(0)
    P_x = Float64(0)
    P_y = Float64(0)   
    
    for j in 1:ny
        for i in 2:nx-1
            if (x1-2 < i && i < x1+rx && y1-1 < j && j < y1+rx)
                # do nothing
            else
                P_x = (p[i,j]-p[i+1,j])/dx[i]
                Dif_x = 0.5/Re*(
                        1/dx[i]*(-(-u[i-1,j]+u[i,j])/dx[i]+(-u[i,j]+u[i+1,j])/dx[i+1])+
                        1/dy[j]*(-(-u[i,jm[j]]+u[i,j])/dya[j]+(-u[i,j]+u[i,jp[j]])/dya[jp[j]]))
                qu[i,j] = P_x + A_x[i,j] + Dif_x + u[i,j]/dt
                sumqu += qu[i,j]^2
            end
            if (x1-1 < i && i < x1+rx && y1-2 < j && j < y1+rx)
                # do nothing
            else
                P_y = (p[i,j]-p[i,jp[j]])/dya[jp[j]]
                Dif_y = 0.5/Re*(
                        1/dx[i]*(-(-v[i-1,j]+v[i,j])/dx[i]+(-v[i,j]+v[i+1,j])/dx[i+1])+
                        1/dya[j]*(-(-v[i,jm[j]]+v[i,j])/dy[j]+(-v[i,j]+v[i,jp[j]])/dy[jp[j]]))
                qv[i,j] = P_y + A_y[i,j] + Dif_y + v[i,j]/dt
                sumqv += qv[i,j]^2

            end 
        end
    end
    return qu,qv,sumqu,sumqv
end

function velocity_revision(nx,ny,rx,dx,dy,dya,qu,qv,sumqu,sumqv,u,v,Re,dt,x1,y1,jp,jm,beta)
    resiu = Float64(0)
    resiv = Float64(0)
    erruv = Float64(0)

    for itr in 1:1000
        sumru = Float64(0)
        sumrv = Float64(0)

        for j in 1:ny
            for i in 2:nx-1             
                if (x1-2 < i && i < x1+rx && y1-1 < j && j < y1+rx)
                    # do nothing
                else
                    resiu = -0.5/Re*(
                        1.0/dx[i]*(-(-u[i-1,j]+u[i,j])/dx[i]+(-u[i,j]+u[i+1,j])/dx[i+1])+
                        1.0/dy[j]*(-(-u[i,jm[j]]+u[i,j])/dya[j]+(-u[i,j]+u[i,jp[j]])/dya[jp[j]])) - qu[i,j] + u[i,j]/dt
                    cx = 1.0/dt + 1.0/dx[i-1]^2 +1.0/dy[j]*(1.0/(2.0*dya[j])+1.0/(2.0*dya[jp[j]]))
                    u[i,j] = u[i,j] - beta*resiu / cx
                    sumru += resiu^2
                end
                if (x1-1 < i && i < x1+rx && y1-2 < j && j < y1+rx)
                    # do nothing
                else
                    resiv = -0.5/Re*(
                        1.0/dx[i]*(-(-v[i-1,j]+v[i,j])/dx[i]+(-v[i,j]+v[i+1,j])/dx[i+1])+
                        1.0/dya[j]*(-(-v[i,jm[j]]+v[i,j])/dy[j]+(-v[i,j]+v[i,jp[j]])/dy[jp[j]])) - qv[i,j] + v[i,j]/dt
                    cy = 1.0/dt + 1.0/dx[i-1]^2 +1.0/dya[j]*(1.0/(2.0*dy[jm[j]])+1.0/(2.0*dy[j]))
                    v[i,j] = v[i,j] - beta*resiv / cy
                    sumrv += resiv^2        
                end
            end
        end
        erruv = sqrt((sumru+sumrv)/(sumqu+sumqv))

        # non-slip condition
        for i in x1:x1+rx-2
            u[i,y1+rx-1] = -u[i,y1+rx]
            u[i,y1]    = -u[i,y1-1]
        end
        for j in y1:y1+rx-2
            v[x1,j] = -v[x1-1,j]
            v[x1+rx-1,j] = -v[x1+rx,j]
        end
        # don't change less than 8
        if(erruv < 1E-8)
            break
        end 
        if(sumqu+sumqv == 0)
            break
        end 
    end
    return u,v,erruv
end

function RHS_Pressure(nx,ny,u,v,dx,dy,dt,x1,y1,rx,jm)
    div = zeros(Float64,nx,ny)
    q = zeros(Float64,nx,ny)
    sumqp = Float64(0)
    
    for j in 1:ny
        for i in 2:nx-1
            if (i > x1-1 && i < x1+rx && j > y1-1 && j < y1+rx)
                # do nothing
            else
                # Consider only fluid part
                div[i,j] = (-u[i-1,j]+u[i,j])/dx[i]+(-v[i,jm[j]]+v[i,j])/dy[j]
                q[i,j] = div[i,j]/dt
                sumqp = sumqp + q[i,j]^2
            end
        end
    end 
    return q,sumqp 
end

function solvepressure(nx,ny,dx,dy,dya,q,sumqp,β,x1,y1,rx,jp,jm)
    errp = Float64(0)
    errpp = Float64(0)
    sumrp = Float64(0)
    resip = Float64(0)
    ss = zeros(Float64,nx,ny)
    
    # 10000 is minimum, but takes much time...
    for itr in 1:10000
        sumrp = 0
        errpp = errp     
        for j in 1:ny                                             
            for i in 2:nx-1
                if (i > x1-1 && i < x1+rx && j > y1-1 && j < y1+rx)
                    # do nothing  
                else
                    resip =   1/(dy[j]*dya[j])*ss[i,jm[j]] +
                              1/dx[i]^2*ss[i-1,j] -
                              (2/dx[i]^2 + 1/(dy[j]*dya[j])+1/(dy[j]*dya[jp[j]]))*ss[i,j] +
                              1/dx[i]^2*ss[i+1,j] +
                              1/(dy[j]*dya[jp[j]])*ss[i,jp[j]] -
                              q[i,j]
                    ss[i,j] = ss[i,j] + β * resip / (2/dx[i]^2 + 1/(dy[j]*dya[j])+1/(dy[j]*dya[jp[j]]))
                    sumrp += resip^2    
                end
            end
        end
        errp = sqrt(sumrp/sumqp)

        # no gradient for inlet
        for j in 1:ny
           ss[ 1,j] = ss[ 2,j]
        end

        # non-slip condition
        for i in x1:x1+rx-1
            ss[i,y1+rx-1] = ss[i,y1+rx]
            ss[i,y1] = ss[i,y1-1]
        end
        for j in y1:y1+rx-1
            ss[x1,j] = ss[x1-1,j]
            ss[x1+rx-1,j] = ss[x1+rx,j]
        end

        if(errp <  1e-6) 
            break
        end 
        if(errpp != 0 && errp*errpp <  1e-12) 
            break
        end 
    end 
    return ss,errp  
end

function Correction(u,v,nx,ny,dt,p,ss,dx,dy,dya,Re,jp,jm)
    for j = 1:ny
        for i = 2:nx-1
            if (x1-2 < i && i < x1+rx && y1-1 < j && j < y1+rx)
            else
                u[i,j] = u[i,j] - dt*(-ss[i,j] + ss[i+1,j])/dx[i]
            end

            if (x1-1 < i && i < x1+rx && y1-2 < j && j < y1+rx)
            else
                v[i,j] = v[i,j] - dt*(-ss[i,j] + ss[i,jp[j]])/dya[jp[j]]
            end
        end
    end

    for j = 1:ny
        for i = 2:nx-1
            p[i,j] = p[i,j] + ss[i,j] - 
            0.5/Re*dt*
            (1/(dy[j]*dya[j])*ss[i,jm[j]] +
            1/dx[i]^2*ss[i-1,j] -
            (2/dx[i]^2 + 1/(dy[j]*dya[j])+1/(dy[j]*dya[jp[j]]))*ss[i,j] +
            1/dx[i]^2*ss[i+1,j] +
            1/(dy[j]*dya[jp[j]])*ss[i,jp[j]])
        end
    end

    # outlet condition
    um = u[nx-1,Int64(ny/2)]  
    for j in 1:ny
        u[nx,j] = u[nx,j] - dt*um*(u[nx,j]-u[nx-1,j])/dx[nx] 
        v[nx,j] = v[nx,j] - dt*um*(v[nx,j]-v[nx-1,j])/dx[nx]  
    end
        
    # non-slip condition
    for i in x1:x1+rx-2
        u[i,y1+rx-1] = -u[i,y1+rx]
        u[i,y1]    = -u[i,y1-1]
    end
    for j in y1:y1+rx-2
        v[x1,j] = -v[x1-1,j]
        v[x1+rx-1,j] = -v[x1+rx,j]
    end

    # set volume-averaged velocity
    for j in 1:ny
        for i in 1:nx
            if i == 1
                up[i,j] = u[i,j]
            elseif (i > x1-1 && i < x1+rx && j > y1-1 && j < y1+rx)
                # do nothing
            else
                up[i,j] = (u[i,j]+u[i-1,j])/2
                vp[i,j] = (v[i,jm[j]]+v[i,j])/2
            end
        end
    end
    return u,v,p,up,vp
end

function q_criterion(nx,ny,jp,jm,dx,dy,dya,u,v)
    Q = zeros(Float64,nx,ny)
    ω = zeros(Float64,nx,ny)
    udx = zeros(Float64,nx,ny)
    vdx = zeros(Float64,nx,ny)
    udy = zeros(Float64,nx,ny)
    vdy = zeros(Float64,nx,ny)

    for j in 1:ny
        for i in 2:nx
            udx[i,j] = (u[i,j]-u[i-1,j])/dx[i]
            vdx[i,j] = (v[i,j]-v[i-1,j])/dx[i]
        end
    end

    for j in 1:ny
        for i in 1:nx
            udy[i,j] = (u[i,j]-u[i,jm[j]])/dya[j]
            vdy[i,j] = (v[i,j]-v[i,jm[j]])/dya[j]
        end
    end

    for j in 1:ny
        for i in 1:nx
            Q[i,j] = -0.5*(udx[i,j]^2+vdy[i,j]^2) - udy[i,j]*vdx[i,j]
            ω[i,j] = vdx[i,j]-udy[i,j]
        end
    end

    return Q,ω
end

function file(nx,ny,up,vp,p,Q,ω,i,xp,yp)
    if i%5 == 0
        open(@sprintf("your_working_dir\\vel_%03d.dat", i),"w") do fileIO
            [(if j1 != ny
                @printf(fileIO, "%11.6G %11.6G %11.6G %11.6G\n",xp[i1],yp[j1],up[i1,j1],vp[i1,j1]) 
            else
                @printf(fileIO, "%11.6G %11.6G %11.6G %11.6G\n\n",xp[i1],yp[j1],up[i1,j1],vp[i1,j1])
            end) for j1=1:ny, i1=1:nx]
        end

        open(@sprintf("your_working_dir\\p_%03d.dat", i),"w") do fileIO
            [(if j1 != ny
                @printf(fileIO, "%11.6G %11.6G %11.6G\n",xp[i1],yp[j1],p[i1,j1]) 
            else
                @printf(fileIO, "%11.6G %11.6G %11.6G\n\n",xp[i1],yp[j1],p[i1,j1])
            end) for j1=1:ny, i1=1:nx]
        end

        open(@sprintf("your_working_dir\\q_%03d.dat", i),"w") do fileIO
            [(if j1 != ny
                @printf(fileIO, "%11.6G %11.6G %11.6G %11.6G\n",xp[i1],yp[j1],Q[i1,j1],ω[i1,j1])
            else
                @printf(fileIO, "%11.6G %11.6G %11.6G %11.6G\n\n",xp[i1],yp[j1],Q[i1,j1],ω[i1,j1])
            end) for j1=1:ny, i1=1:nx]
        end

    end
end

### main program ###
begin
    # simulation parameters
    const Re = Float64(100)
    const dt = Float64(5e-3)
    const time_iteration = Int64(10000)

    # mesh parameters 
    const nx = Int64(768)
    const ny = Int64(256)
    #const ALG = Float64(0.95)           # only for un-uniform mesh
    
    # Rectangle information
    const x1 = Int64(128)                # Rectangle start point x
    const rx = Int64(20)                 # Rectangle size    
    const y1 = Int64((ny-rx)/2+1)        # Rectangle start point y

    xp = collect(Float64, 0:0.05:38.4)
    dx = fill(Float64(0.05), nx)
    yp = collect(Float64, 0:0.05:12.8)
    dy = fill(Float64(0.05), ny)
    dya = fill(Float64(0.05), ny) 

    jp,jm = periodic(ny)
            
    # initialisation
    u = ones(Float64,nx,ny)
    for j = y1-1:y1+rx
        for i = x1-1:x1+rx
            u[i,j] = 0
        end
    end
    v = zeros(Float64,nx,ny)
    p = zeros(Float64,nx,ny) 
    Q = zeros(Float64,nx,ny)
    ω = zeros(Float64,nx,ny)           
    up = zeros(Float64,nx,ny)
    vp = zeros(Float64,nx,ny)

    Con_x_p = zeros(Float64,nx,ny)
    Con_y_p = zeros(Float64,nx,ny)

   ##### Begin time-stepping #####
    println("start")
    i = Int64(0)
    t = Float64(0)

    while i < time_iteration
        global t += dt
        global i += 1
        
        global A_x, A_y, Con_x_p, Con_y_p  = convection(nx,ny,rx,dx,dy,dya,u,v,Con_x_p,Con_y_p,jp,jm)
        global qu,qv,sumqu,sumqv = diffusion(nx,ny,rx,dx,dy,dya,Re,dt,u,v,A_x, A_y,p,x1,y1,jp,jm)
        global u,v,erruv = velocity_revision(nx,ny,rx,dx,dy,dya,qu,qv,sumqu,sumqv,u,v,Re,dt,x1,y1,jp,jm,1.6)
        global q,sumqp = RHS_Pressure(nx,ny,u,v,dx,dy,dt,x1,y1,rx,jm)
        global ss,errp = solvepressure(nx,ny,dx,dy,dya,q,sumqp,1.6,x1,y1,rx,jp,jm)
        global u,v,p,up,vp = Correction(u,v,nx,ny,dt,p,ss,dx,dy,dya,Re,jp,jm)
        global Q,ω = q_criterion(nx,ny,jp,jm,dx,dy,dya,u,v)
        file(nx,ny,up,vp,p,Q,ω,i,xp,yp)
       
        @printf("Step:%05d, time:%0.3G, erruv:%10.8G, errp: %10.8G\n",i,t,erruv,errp)

        # in case the result is apparently wrong, this program will be finished
        if(u[10,10]>1e2)
            println("break!")
            break
        end
    end
end
