tic
heat();
results1(ii) = toc;

tic
heat_Vec();
results2(ii) = toc;

tic
heat_wof();
results3(ii) = toc;

tic
heat_Vec_wof();
results4(ii) = toc;

function heat()
    N = 251;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = zeros(N);
    y = zeros(N);
    for j=1:N
        for i=1:N
            x(i,j) = (i-1)*dx;
        end
    end
    for j=1:N
        for i=1:N
            y(i,j) = (j-1)*dx;
        end
    end
    dt = dx;
    tmax = (10-mod(10,dt))/dt;
    k = 0.25*dx;
    alpha = k*dt/dx^2;
    T0 = zeros(N);
    T1 = zeros(N);
    for i=1:N
        T1(1,i) = 1;
        T1(end,i) = 1;
        T1(i,1) = 1;
        T1(i,end) = 1;
    end

    for i=1:N
        T0(1,i) = 1;
        T0(end,i) = 1;
        T0(i,1) = 1;
        T0(i,end) = 1;
    end
    t = 0*dt;
    fileID = fopen('MATLAB.dat','wt');
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "y", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
    for j=1:N
        for i=1:N
            fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
        end
    end

    for k=1:tmax
        for j= 2:N-1
            for i=2:N-1
                T1(i,j)=T0(i,j)+alpha*((T0(i+1,j)-2*T0(i,j)+T0(i-1,j))+(T0(i,j+1)-2*T0(i,j)+T0(i,j-1)));
            end
        end
        T0=T1;
        if mod(k,500)==0
            t = k*dt;
            fprintf(fileID,'title ="ZoneTime_%f"\n',t);
            fprintf(fileID,'variables = "x", "y", "u",\n');
            fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
            for j=1:N
                for i=1:N
                    fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
                end
            end
        end
    end
    t = tmax*dt;
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "y", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
    for j=1:N
        for i=1:N
            fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
        end
    end
    fclose(fileID);
end

function heat_Vec()
    N = 251;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    x = zeros(N);
    y = zeros(N);
    x(:,1) = (0:N-1)*dx;
    y(1,:) = (0:N-1)*dx;
    for i=1:N
        x(:,i) = x(:,1);
        y(i,:) = y(1,:);
    end

    dt = dx;
    tmax = (10-mod(10,dt))/dt;
    k = 0.25*dx;
    alpha = k*dt/dx^2;
    T0 = zeros(N);
    T1 = zeros(N);
    T1(1,:) = 1;
    T1(end,:) = 1;
    T1(:,1) = 1;
    T1(:,end) = 1;

    T0(1,:) = 1;
    T0(end,:) = 1;
    T0(:,1) = 1;
    T0(:,end) = 1;
    t = 0*dt;
    fileID = fopen('MATLAB.dat','wt');
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "y", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
    for j=1:N
        for i=1:N
            fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
        end
    end

    for k=1:tmax
        T1(2:N-1,2:N-1)=T0(2:N-1,2:N-1)+alpha*((T0(3:N,2:N-1)-2*T0(2:N-1,2:N-1)+T0(1:N-2,2:N-1))+...
                                               (T0(2:N-1,3:N)-2*T0(2:N-1,2:N-1)+T0(2:N-1,1:N-2)));
        T0=T1;
        if mod(k,500)==0
            t = k*dt;
            fprintf(fileID,'title ="ZoneTime_%f"\n',t);
            fprintf(fileID,'variables = "x", "y", "u",\n');
            fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
            for j=1:N
                for i=1:N
                    fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
                end
            end
        end
    end
    t = tmax*dt;
    fprintf(fileID,'title ="ZoneTime_%f"\n',t);
    fprintf(fileID,'variables = "x", "y", "u",\n');
    fprintf(fileID,'zone T="Zone_%f" i=%d j=%d\n',t,N,N);
    for j=1:N
        for i=1:N
            fprintf(fileID,'%1.9e %1.9e %1.9e\n',x(i,j),y(i,j),T0(i,j));
        end
    end
    fclose(fileID);
end

function heat_wof()
    N = 251;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);
    dt = dx;
    tmax = (10-mod(10,dt))/dt;
    k = 0.25*dx;
    alpha = k*dt/dx^2;
    T0 = zeros(N);
    T1 = zeros(N);
    for i=1:N
        T1(1,i) = 1;
        T1(end,i) = 1;
        T1(i,1) = 1;
        T1(i,end) = 1;
    end

    for i=1:N
        T0(1,i) = 1;
        T0(end,i) = 1;
        T0(i,1) = 1;
        T0(i,end) = 1;
    end
    
    for k=1:tmax
        for j= 2:N-1
            for i=2:N-1
                T1(i,j)=T0(i,j)+alpha*((T0(i+1,j)-2*T0(i,j)+T0(i-1,j))+(T0(i,j+1)-2*T0(i,j)+T0(i,j-1)));
            end
        end
        T0=T1;
        
    end
    
end

function heat_Vec_wof()
    N = 251;
    a = 0;
    b = pi;
    dx = (b-a)/(N-1);

    dt = dx;
    tmax = (10-mod(10,dt))/dt;
    k = 0.25*dx;
    alpha = k*dt/dx^2;
    T0 = zeros(N);
    T1 = zeros(N);
    
    T1(1,:) = 1;
    T1(end,:) = 1;
    T1(:,1) = 1;
    T1(:,end) = 1;

    T0(1,:) = 1;
    T0(end,:) = 1;
    T0(:,1) = 1;
    T0(:,end) = 1;

    for k=1:tmax
        T1(2:N-1,2:N-1)=T0(2:N-1,2:N-1)+alpha*((T0(3:N,2:N-1)-2*T0(2:N-1,2:N-1)+T0(1:N-2,2:N-1))+...
                                               (T0(2:N-1,3:N)-2*T0(2:N-1,2:N-1)+T0(2:N-1,1:N-2)));
        T0=T1;
    end
end