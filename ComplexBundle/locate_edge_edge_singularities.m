function [X_EE1,X_EE2,Inter_index,E2V_1,E2V_2] = locate_edge_edge_singularities(z, M1, M2, omega_EV, omega_VE, sing_EE, tol)
idee = find(abs(sing_EE(:)) > 0.5);

nout = 0;
Inter_index = cast(zeros(M1.ne,M2.ne), 'uint32');
X_EE1 = zeros(length(idee),3);
X_EE2 = zeros(length(idee),3);
X2E = zeros(length(idee),2);
X2bar = zeros(length(idee),2);
for i = 1:length(idee)
    [ide1,ide2] = ind2sub([M1.ne,M2.ne], idee(i));
    idvx1 = M1.E2V(ide1,:);
    idvx2 = M2.E2V(ide2,:);
    li = abs(z(idvx1,idvx2));
    li = li/mean(li(:));
    om1 = omega_EV(ide1,idvx2);
    om2 = omega_VE(idvx1,ide2);
    assert(abs(abs(om1(1) - om1(2) + om2(2) - om2(1)) - 2*pi) < tol);

    a_bl = li(1,1) - li(2,1)*exp(1i*om1(1)) - li(1,2)*exp(1i*om2(1)) + li(2,2)*exp(1i*om1(2) + 1i*om2(1));
    b_bl =-li(1,1) + li(1,2)*exp(1i*om2(1));
    c_bl =-li(1,1) + li(2,1)*exp(1i*om1(1));
    d_bl = li(1,1);
    Fct = @(X,Y) X.*Y*a_bl + X*b_bl + Y*c_bl + d_bl;
    
    a = imag(b_bl)*real(a_bl) - real(b_bl)*imag(a_bl);
    b = imag(d_bl)*real(a_bl) - real(d_bl)*imag(a_bl) + imag(b_bl)*real(c_bl) - real(b_bl)*imag(c_bl);
    c = imag(d_bl)*real(c_bl) - real(d_bl)*imag(c_bl);
    get_y_from_x = @(x)-(abs(real(a_bl*x + c_bl)) >= abs(imag(a_bl*x + c_bl)))*real(b_bl*x + d_bl)/real(a_bl*x + c_bl) + ...
                       -(abs(real(a_bl*x + c_bl)) <  abs(imag(a_bl*x + c_bl)))*imag(b_bl*x + d_bl)/imag(a_bl*x + c_bl);

    if abs(a) > tol
        Delta = b^2 - 4*a*c;
        if Delta < 0
            error('No singularity');
        end

        x1 = (- b + sqrt(Delta))/(2*a);
        y1 = get_y_from_x(x1);
        x2 = (- b - sqrt(Delta))/(2*a);
        y2 = get_y_from_x(x2);
        if x1 >= 0 && x1 <= 1 && y1 >= 0 && y1 <= 1
            le = [x1; y1];
        elseif x2 >= 0 && x2 <= 1 && y2 >= 0 && y2 <= 1
            le = [x2; y2];
        else
            warning('Singularity outside of edge-edge face.');
        end
    else
        le(1) =-c/b;
        le(2) = get_y_from_x(le(1));
    end

    if any(le < 0) || any(le > 1) || abs(Fct(le(1), le(2))) > tol
        nout = nout + 1;

        [X,Y] = meshgrid(0:0.01:1, 0:0.01:1);
        F =  Fct(X, Y);
        figure; quiver(X, Y, real(F), imag(F));
        hold on; plot(le(1), le(2), 'rx'); hold off; axis equal;
        pause
    end
    assert(abs(Fct(le(1), le(2))) < tol, ['Could not find singularity. Err: ', num2str(abs(Fct(le(1), le(2)))), ' -- tol: ', num2str(tol)]);

    Inter_index(ide1,ide2) = M1.nv + M2.nv + i;
    X_EE1(i,:) = M1.X(idvx1(1),:)*(1-le(2)) + M1.X(idvx1(2),:)*le(2);
    X_EE2(i,:) = M2.X(idvx2(1),:)*(1-le(1)) + M2.X(idvx2(2),:)*le(1);
    X2E(i,:) = [ide1, ide2];
    X2bar(i,:) = [le(2), le(1)];
end
disp(['Number of singu outside edge-edge: ', num2str(sum(nout))]);

% Track edges from initial meshes
E2V_1 = cell(M1.ne,1);
for i = 1:M1.ne
    idvx = find(X2E(:,1) == i);
    [~,id] = sort(X2bar(idvx,1));
    idvx = [M1.E2V(i,1); M1.nv + M2.nv + idvx(id); M1.E2V(i,2)];
    E2V_1{i} = [i*ones(length(idvx)-1,1), idvx(1:end-1), idvx(2:end)];
end
E2V_1 = cell2mat(E2V_1);

E2V_2 = cell(M2.ne,1);
for i = 1:M2.ne
    idvx = find(X2E(:,2) == i);
    [~,id] = sort(X2bar(idvx,2));
    idvx = [M1.nv + M2.E2V(i,1); M1.nv + M2.nv + idvx(id); M1.nv + M2.E2V(i,2)];
    E2V_2{i} = [i*ones(length(idvx)-1,1), idvx(1:end-1), idvx(2:end)];
end
E2V_2 = cell2mat(E2V_2);
end