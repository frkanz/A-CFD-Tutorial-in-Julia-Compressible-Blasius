using Printf
using BenchmarkTools
function heat()
    N = 1001
    T0 = Matrix{Float64}(undef,N,N)
    T1 = Matrix{Float64}(undef,N,N)
    x  = Matrix{Float64}(undef,N,N)
    y  = Matrix{Float64}(undef,N,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    for j=1:N
        for i=1:N
            x[i,j] = (i-1)*dx
        end
    end
    for j=1:N
        for i=1:N
            y[i,j] = (j-1)*dx
        end
    end

    dt = dx
    tmax = (10-mod(10,dt))/dt
    k = 0.25*dx
    alpha = k*dt/dx^2
    for i=1:N
        T1[1,i] = 1.0
        T1[N,i] = 1.0
        T1[i,1] = 1.0
        T1[i,N] = 1.0
    end
    for i=1:N
        T0[1,i] = 1.0
        T0[N,i] = 1.0
        T0[i,1] = 1.0
        T0[i,N] = 1.0
    end
    t = 0*dt
    vel = open("Julia.dat", "w")
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"y\", \"T\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
    for j=1:N
        for i=1:N
            write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
        end
    end

    for k=1:tmax
        for j=2:N-1
            for i=2:N-1
                T1[i,j] = T0[i,j]+alpha*((T0[i+1,j]-2*T0[i,j]+T0[i-1,j])+(T0[i,j+1]-2*T0[i,j]+T0[i,j-1]))
            end
        end
        T0 .= T1
        if mod(k,500)==0
            write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
            write(vel, "variables = \"x\", \"y\", \"T\""," \n")
            write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
            for j=1:N
                for i=1:N
                    write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
                end
            end
        end
    end
    t = tmax*dt
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"y\", \"T\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
    for j=1:N
        for i=1:N
            write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
        end
    end
    close(vel)
end

function heat_Vec()
    N = 1001
    T0 = Matrix{Float64}(undef,N,N)
    T1 = Matrix{Float64}(undef,N,N)
    x  = Matrix{Float64}(undef,N,N)
    y  = Matrix{Float64}(undef,N,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x[:,1] = (0:N-1)*dx
    x[:,:] .= x[:,1]
    y[1,:] = (0:N-1)*dx
    y[:,:] .= y[1,:]

    dt = dx
    tmax = (10-mod(10,dt))/dt
    k = 0.25*dx
    alpha = k*dt/dx^2
    T1[1,:] .= 1.0
    T1[N,:] .= 1.0
    T1[:,1] .= 1.0
    T1[:,N] .= 1.0

    T0[1,:] .= 1.0
    T0[N,:] .= 1.0
    T0[:,1] .= 1.0
    T0[:,N] .= 1.0
    t = 0*dt
    vel = open("Julia.dat", "w")
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"y\", \"T\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
    for j=1:N
        for i=1:N
            write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
        end
    end

    for k=1:tmax
        T1[2:N-1,2:N-1] .= T0[2:N-1,2:N-1].+alpha.*((T0[3:N,2:N-1].-2 .*T0[2:N-1,2:N-1].+T0[1:N-2,2:N-1]).+
                                                    (T0[2:N-1,3:N].-2 .*T0[2:N-1,2:N-1].+T0[2:N-1,1:N-2]))
        T0 .= T1
        if mod(k,500)==0
            write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
            write(vel, "variables = \"x\", \"y\", \"T\""," \n")
            write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
            for j=1:N
                for i=1:N
                    write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
                end
            end
        end
    end
    t = tmax*dt
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"y\", \"T\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," j=",@sprintf("%d",N)," \n")
    for j=1:N
        for i=1:N
            write(vel, @sprintf("%1.9e",x[i,j])," ",@sprintf("%1.9e",y[i,j])," ",@sprintf("%1.9e",T0[i,j])," \n")
        end
    end
    close(vel)
end

function heat_wof()
    N = 1001
    T0 = Matrix{Float64}(undef,N,N)
    T1 = Matrix{Float64}(undef,N,N)
    x  = Matrix{Float64}(undef,N,N)
    y  = Matrix{Float64}(undef,N,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    for j=1:N
        for i=1:N
            x[i,j] = (i-1)*dx
        end
    end
    for j=1:N
        for i=1:N
            y[i,j] = (j-1)*dx
        end
    end

    dt = dx
    tmax = (10-mod(10,dt))/dt
    k = 0.25*dx
    alpha = k*dt/dx^2
    for i=1:N
        T1[1,i] = 1.0
        T1[N,i] = 1.0
        T1[i,1] = 1.0
        T1[i,N] = 1.0
    end
    for i=1:N
        T0[1,i] = 1.0
        T0[N,i] = 1.0
        T0[i,1] = 1.0
        T0[i,N] = 1.0
    end

    for k=1:tmax
        for j=2:N-1
            for i=2:N-1
                T1[i,j] = T0[i,j]+alpha*((T0[i+1,j]-2*T0[i,j]+T0[i-1,j])+(T0[i,j+1]-2*T0[i,j]+T0[i,j-1]))
            end
        end
        T0 .= T1
    end
end

function heat_Vec_wof()
    N = 1001
    T0 = Matrix{Float64}(undef,N,N)
    T1 = Matrix{Float64}(undef,N,N)
    x  = Matrix{Float64}(undef,N,N)
    y  = Matrix{Float64}(undef,N,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x[:,1] = (0:N-1)*dx
    x[:,:] .= x[:,1]
    y[1,:] = (0:N-1)*dx
    y[:,:] .= y[1,:]

    dt = dx
    tmax = (10-mod(10,dt))/dt
    k = 0.25*dx
    alpha = k*dt/dx^2
    T1[1,:] .= 1.0
    T1[N,:] .= 1.0
    T1[:,1] .= 1.0
    T1[:,N] .= 1.0

    T0[1,:] .= 1.0
    T0[N,:] .= 1.0
    T0[:,1] .= 1.0
    T0[:,N] .= 1.0

    for k=1:tmax
        T1[2:N-1,2:N-1] .= T0[2:N-1,2:N-1].+alpha.*((T0[3:N,2:N-1].-2 .*T0[2:N-1,2:N-1].+T0[1:N-2,2:N-1]).+
                                                    (T0[2:N-1,3:N].-2 .*T0[2:N-1,2:N-1].+T0[2:N-1,1:N-2]))
        T0 .= T1
    end

end

results1 = @benchmark heat()
results2 = @benchmark heat_Vec()
results3 = @benchmark heat_wof()
results4 = @benchmark heat_Vec_wof()
