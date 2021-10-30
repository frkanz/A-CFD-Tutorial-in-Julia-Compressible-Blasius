tic
burgers();
results1 = toc;

tic
burgers_Vec();
results2 = toc;

tic
burgers_wof();
results3 = toc;

tic
burgers_Vec_wof();
results4 = toc;

function burgers()
    N = 2500;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = a:dx:b;
    dt = 0.5*dx;
    tmax = (1-mod(1,dt))/dt;
    U = zeros(N,1);
    u = zeros(N,1);
    for i=1:N
        u(i) = sin(x(i));
    end
    fileID = fopen('MATLAB.dat','wt');
    t = 0*dt;
    U(1) = u(1);
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
    for i=1:N
        fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
    end

    for j=1:tmax
        for i=2:N-1
            U(i) = u(i)-dt/dx*(0.5*u(i)^2-0.5*(u(i-1)^2));
        end
        u = U;
        if mod(j,100)==0
            t = j*dt;
            fprintf(fileID,'title ="ZoneTime_%f"\n',t);
            fprintf(fileID,'variables = "x", "u",\n');
            fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
            for i=1:N
                fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
            end
        end
    end
    t = tmax*dt;
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
    for i=1:N
        fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
    end
    fclose(fileID);
end

function burgers_Vec()
    N = 2500;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = a:dx:b;
    dt = 0.5*dx;
    tmax = (1-mod(1,dt))/dt;
    U = zeros(N,1);
    u = zeros(N,1);
    u = sin(x);
    fileID = fopen('MATLAB.dat','wt');
    t = 0*dt;
    U(1) = u(1);
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
    for i=1:N
        fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
    end

    for j=1:tmax
        U(2:N-1) = u(2:N-1)-dt/dx.*(0.5.*u(2:N-1).^2-0.5.*(u(1:N-2).^2));
        u = U;
        if mod(j,100)==0
            t = j*dt;
            fprintf(fileID,'title ="ZoneTime_%f"\n',t);
            fprintf(fileID,'variables = "x", "u",\n');
            fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
            for i=1:N
                fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
            end
        end
    end
    t = tmax*dt;
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d\n',t,N+1);
    for i=1:N
        fprintf(fileID,'%1.9e %1.9e\n',x(i),u(i));
    end
    fclose(fileID);
end

function burgers_wof()
    N = 2500;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = a:dx:b;
    dt = 0.5*dx;
    tmax = (1-mod(1,dt))/dt;
    U = zeros(N,1);
    u = zeros(N,1);
    for i=1:N
        u(i) = sin(x(i));
    end
    U(1) = u(1);
    for j=1:tmax
        for i=2:N-1
            U(i) = u(i)-dt/dx*(0.5*u(i)^2-0.5*(u(i-1)^2));
        end
        u = U;
    end
end

function burgers_Vec_wof()
    N = 2500;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = a:dx:b;
    dt = 0.5*dx;
    tmax = (1-mod(1,dt))/dt;
    U = zeros(N,1);
    u = zeros(N,1);
    u = sin(x);
    U(1) = u(1);
    for j=1:tmax
        U(2:N-1) = u(2:N-1)-dt/dx.*(0.5.*u(2:N-1).^2-0.5.*(u(1:N-2).^2));
        u = U;
    end
end
