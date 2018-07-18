subplot(2,3,1)
tri = delaunay(x(:,1),x(:,2));
trisurf(tri, x(:,1),x(:,2),y)
xlabel('theta1')
ylabel('theta2')

subplot(2,3,2)
tri = delaunay(x(:,1),x(:,3));
trisurf(tri, x(:,1),x(:,3),y)
xlabel('theta1')
ylabel('theta3')

subplot(2,3,3)
tri = delaunay(x(:,1),x(:,4));
trisurf(tri, x(:,1),x(:,4),y)
xlabel('theta1')
ylabel('theta4')

subplot(2,3,4)
tri = delaunay(x(:,2),x(:,3));
trisurf(tri, x(:,2),x(:,3),y)
xlabel('theta2')
ylabel('theta3')

subplot(2,3,5)
tri = delaunay(x(:,2),x(:,4));
trisurf(tri, x(:,2),x(:,4),y)
xlabel('theta2')
ylabel('theta4')

subplot(2,3,6)
tri = delaunay(x(:,3),x(:,4));
trisurf(tri, x(:,3),x(:,4),y)
xlabel('theta3')
ylabel('theta4')