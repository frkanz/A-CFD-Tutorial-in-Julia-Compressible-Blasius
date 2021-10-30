using Printf
using BenchmarkTools
function burgers()
    N = 10001
    U = Vector{Float64}(undef,N)
    u = Vector{Float64}(undef,N)
    x = Vector{Float64}(undef,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x .= a:dx:b
    dt = 0.5*dx
    tmax = (1-mod(1,dt))/dt
    for i=1:N
        u[i] = sin(x[i])
    end
    U[1] = u[1]
    t = 0*dt
    vel = open("Julia.dat", "w")
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"u\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
    for i=1:N
        write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
    end
    for j=1:tmax
        for i=2:N-1
            U[i] = u[i]-dt/dx*(0.5*u[i]^2-0.5*(u[i-1]^2))
        end
        u .= U
        if mod(j,100)==0
            t = j*dt
            write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
            write(vel, "variables = \"x\", \"u\""," \n")
            write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
            for i=1:N
                write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
            end
        end
    end
    t = tmax*dt
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"u\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
    for i=1:N
        write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
    end
    close(vel)
end

function burgers_Vec()
    N = 10001
    U = Vector{Float64}(undef,N)
    u = Vector{Float64}(undef,N)
    x = Vector{Float64}(undef,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x .= a:dx:b
    dt = 0.5*dx
    tmax = (1-mod(1,dt))/dt
    u .= sin.(x)
    U[1] = u[1]
    t = 0*dt
    vel = open("Julia.dat", "w")
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"u\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
    for i=1:N
        write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
    end
    for j=1:tmax
        U[2:N-1] .= u[2:N-1].-dt./dx.*(0.5.*u[2:N-1].^2 .- 0.5.*(u[1:N-2].^2))
        u .= U
        if mod(j,100)==0
            write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
            write(vel, "variables = \"x\", \"u\""," \n")
            write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
            for i=1:N
                write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
            end
        end
    end
    t = tmax*dt
    write(vel, "title =\"ZoneTime_",string(t),"\""," \n")
    write(vel, "variables = \"x\", \"u\""," \n")
    write(vel, "zone T=\"Zone_"*string(t)*"\" i=",@sprintf("%d",N)," \n")
    for i=1:N
        write(vel, @sprintf("%1.9e",x[i])," ",@sprintf("%1.9e",u[i])," \n")
    end
    close(vel)
end

function burgers_wof()
    N = 10001
    U = Vector{Float64}(undef,N)
    u = Vector{Float64}(undef,N)
    x = Vector{Float64}(undef,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x .= a:dx:b
    dt = 0.5*dx
    tmax = (1-mod(1,dt))/dt
    for i=1:N
        u[i] = sin(x[i])
    end
    U[1] = u[1]
    for j=1:tmax
        for i=2:N-1
            U[i] = u[i]-dt/dx*(0.5*u[i]^2-0.5*u[i-1]^2)
        end
        u .= U
    end
end

function burgers_Vec_wof()
    N = 10001
    U = Vector{Float64}(undef,N)
    u = Vector{Float64}(undef,N)
    x = Vector{Float64}(undef,N)
    a = 0
    b = π;
    dx = (b-a)/(N-1)
    x .= a:dx:b
    dt = 0.5*dx
    tmax = (1-mod(1,dt))/dt
    u .= sin.(x)
    U[1] = u[1]
    for j=1:tmax
        U[2:N-1] .= u[2:N-1].-dt./dx.*(0.5.*u[2:N-1].^2 .- 0.5.*(u[1:N-2].^2))
        u .= U
    end
end

results1 = @benchmark burgers()
results2 = @benchmark burgers_Vec()
results3 = @benchmark burgers_wof()
results4 = @benchmark burgers_Vec_wof()
