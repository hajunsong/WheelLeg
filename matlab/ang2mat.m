function mat = ang2mat(psi, theta, phi)

    Rz1 = [cos(psi), -sin(psi), 0;
        sin(psi), cos(psi), 0;
        0, 0, 1];
    Rx = [1, 0, 0;
        0, cos(theta), -sin(theta);
        0, sin(theta), cos(theta)];
    Rz2 = [cos(phi), -sin(phi), 0;
        sin(phi), cos(phi), 0;
        0, 0, 1];

    mat = Rz1*Rx*Rz2;

end