struct Sphere{
    vec3 c;// Center
    float r;// Radius
    int i;// Texture Id
};

struct Ellipsoid {
    vec3 c;  // Center
    vec3 radii;  // Radii along each axis (x, y, z)
    int i;  // Texture Id
};

struct Plane{
    vec3 n;// Normal
    vec3 p;// Point
    int i;// Texture Id
};

struct Hit{
    float t;// Intersection depth
    vec3 n;// Normal
    int i;// Texture Id
};

struct Ray{
    vec3 o;// Origin
    vec3 d;// Direction
};

struct Material
{
    vec3 d;// Diffuse
};

struct Box{
    vec3 mini;
    vec3 maxi;
    int i;
};

struct Cylinder{
    vec3 c; 
    float r;
    float h;
    int i;
};

struct Capsule {
    vec3 start;  // Point de départ d'une extrémité de la pilule
    vec3 end;    // Point d'arrêt de l'autre extrémité de la pilule
    float radius; // Rayon de la pilule
    int i;        // Identifiant de texture de la pilule
};



float Checkers(in vec2 p)
{
    // Filter kernel
    vec2 w=fwidth(p)+.001;
    // Box box filter
    vec2 i=2.*(abs(fract((p-.5*w)*.5)-.5)-abs(fract((p+.5*w)*.5)-.5))/w;
    // xor pattern
    return.5-.5*i.x*i.y;
}

// Compute point on ray
vec3 Point(Ray ray,float t)
{
    return ray.o+t*ray.d;
}

// Compute color
// i : Texture index
// p : Point
Material Texture(vec3 p,int i)
{
    if(i==1)
    {
        return Material(vec3(.8,.5,.4));
    }
    else if(i==0)
    {
        // compute checkboard
        float f=Checkers(.5*p.xy);
        vec3 col=vec3(.4,.5,.7)+f*vec3(.1);
        return Material(col);
    }
    return Material(vec3(0));
}

// Sphere intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectSphere(Ray ray,Sphere sph,out Hit x)
{
    vec3 oc=ray.o-sph.c;
    float b=dot(oc,ray.d);
    float c=dot(oc,oc)-sph.r*sph.r;
    float d=b*b-c;
    if(d>0.)
    {
        float t=-b-sqrt(d);
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            x=Hit(t,normalize(p-sph.c),sph.i);
            
            return true;
        }
    }
    return false;
    
}

// Plane intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectPlane(Ray ray,Plane pl,out Hit x)
{
    float t=-dot(ray.o-pl.p,pl.n)/dot(ray.d,pl.n);
    if(t>0.)
    {
        
        x=Hit(t,vec3(0,0,1),0);
        return true;
    }
    return false;
}

bool IntersectBox(Ray ray, Box box, out Hit x){

    float tMin = -99999.0; //distance d entree la plus proche
    float tMax = 99999.0; //distance de sortie la plus lointaine
    
    //on calcule les potentielles distances d'entrées et de sortie
    //pour chaque dimensions
    for(int i = 0; i < 3; i++){
    
        float t0 = (box.mini[i] - ray.o[i])/ray.d[i];
        float t1 = (box.maxi[i] - ray.o[i])/ray.d[i];
        
        if(t1 < t0){
           
            float tmp = t1;
            t1 = t0;
            t0 = tmp;
        }
        
        tMin = max(tMin,t0);
        tMax = min(tMax, t1);
        
        if(tMin > tMax){
            
            return false;
        }
    }
    
    vec3 p=Point(ray,tMin);
    x=Hit(tMin,normalize(p-(box.mini + box.maxi) * 0.5),box.i);
    


    return true;
}

bool IntersectEllipsoid(Ray ray, Ellipsoid ellipsoid, out Hit x) {
    // Transform the ray into ellipsoid space (inverse scale and translation)
    vec3 oc = (ray.o - ellipsoid.c) / ellipsoid.radii;
    vec3 rd = ray.d / ellipsoid.radii;

    float a = dot(rd, rd);
    float b = dot(oc, rd);
    float c = dot(oc, oc) - 1.0;

    float discriminant = b * b - a * c;

    if (discriminant > 0.0) {
        float t1 = (-b - sqrt(discriminant)) / a;
        float t2 = (-b + sqrt(discriminant)) / a;

        if (t1 > 0.0 || t2 > 0.0) {
            float t = (t1 > 0.0) ? t1 : t2;
            vec3 p = Point(ray, t);
            vec3 normal = normalize(vec3(2.0 * p.x, 2.0 * p.y, 2.0 * p.z));
            x = Hit(t, normal, ellipsoid.i);
            return true;
        }
    }

    return false;
}

bool IntersectCylinder(Ray ray, Cylinder cy, out Hit h){
    vec3 rayOriginLocal = ray.o - cy.c;

    float a = ray.d.x*ray.d.x + ray.d.y*ray.d.y;
    float b = 2.0 * (ray.d.x * rayOriginLocal.x + ray.d.y * rayOriginLocal.y);
    float c = rayOriginLocal.x * rayOriginLocal.x + rayOriginLocal.y * rayOriginLocal.y - cy.r * cy.r;

    float delta = b * b - 4.0 * a * c;

    if (delta > 0.0){
        float t1 = (-b - sqrt(delta)) / (2.0 * a);
        float t2 = (-b + sqrt(delta)) / (2.0 * a);
        float t;
        
        if (t1 < t2) t = t1;
        else t = t2;
        
        float intersectXLocal = ray.d.x * t + rayOriginLocal.x;
        float intersectYLocal = ray.d.y * t + rayOriginLocal.y;
        float intersectZLocal = ray.d.z * t + rayOriginLocal.z;
        
        if (intersectXLocal >= -cy.h / 2.0 && intersectXLocal <= cy.h / 2.0){
                                
                vec3 p=Point(ray,t);
                h=Hit(t,normalize(p-cy.c),cy.i);

                return true;
            }
    }

    return false;   
}

// Fonction pour résoudre une équation quadratique ax^2 + bx + c = 0
bool SolveQuadratic(float a, float b, float c, out float x1, out float x2) {
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant > 0.0) {
        float sqrtDiscriminant = sqrt(discriminant);
        x1 = (-b - sqrtDiscriminant) / (2.0 * a);
        x2 = (-b + sqrtDiscriminant) / (2.0 * a);
        return true;  // Deux solutions réelles
    } else if (discriminant == 0.0) {
        x1 = -b / (2.0 * a);
        x2 = x1;
        return true;  // Une seule solution réelle
    } else {
        x1 = 0.0;
        x2 = 0.0;
        return false;  // Pas de solution réelle
    }
}


bool IntersectCapsule(Ray ray, Capsule capsule, out Hit x) {
    // Intersection avec les sphères aux extrémités
    Sphere sphere1 = Sphere(capsule.start, capsule.radius, capsule.i);
    Sphere sphere2 = Sphere(capsule.end, capsule.radius, capsule.i);
    Hit hit1, hit2;
    bool intersects1 = IntersectSphere(ray, sphere1, hit1);
    bool intersects2 = IntersectSphere(ray, sphere2, hit2);

    // Intersection avec le cylindre
    vec3 dir = capsule.end - capsule.start;
    float cylinderLength = length(dir);
    dir = normalize(dir);
    vec3 oc1 = ray.o - capsule.start;
    vec3 oc2 = ray.o - capsule.end;

    float b1 = dot(oc1, dir);
    float b2 = dot(oc2, dir);
    float c1 = dot(oc1, oc1) - capsule.radius * capsule.radius;
    float c2 = dot(oc2, oc2) - capsule.radius * capsule.radius;

    float t1, t2;
    bool intersectsCylinder1 = SolveQuadratic(1.0, 2.0 * b1, c1, t1, t2);
    bool intersectsCylinder2 = SolveQuadratic(1.0, 2.0 * b2, c2, t1, t2);

    // Trouver la plus proche intersection
    bool intersectsCylinder = false;
    float t = 0.0;

    if (intersectsCylinder1 && (t1 > 0.0 && t1 < cylinderLength)) {
        intersectsCylinder = true;
        t = t1;
    }
    if (intersectsCylinder2 && (t2 > 0.0 && t2 < cylinderLength)) {
        if (!intersectsCylinder || t2 < t) {
            intersectsCylinder = true;
            t = t2;
        }
    }

    // Trouver la plus proche intersection parmi toutes les parties de la capsule
    bool intersects = false;
    x = Hit(1000.0, vec3(0), -1);

    if (intersects1 && (!intersects || hit1.t < x.t)) {
        intersects = true;
        x = hit1;
    }
    if (intersects2 && (!intersects || hit2.t < x.t)) {
        intersects = true;
        x = hit2;
    }
    if (intersectsCylinder && (!intersects || t < x.t)) {
        intersects = true;
        x = Hit(t, normalize(Point(ray, t) - (capsule.start + capsule.end) * 0.5), capsule.i);
    }

    return intersects;
}




// Scene intersection
// ray : The ray
//   x : Returned intersection information
bool Intersect(Ray ray,out Hit x)
{
    // Spheres
    //const Sphere sph1=Sphere(vec3(3.5,0.,3.),1.,1);
    const Sphere sph2=Sphere(vec3(3,0.,6.),1.,1);
    
    // Plane
    const Plane pl=Plane(vec3(0.,0.,1.),vec3(0.,0.,0.),0);
    
    // Box
    const Box bx=Box(vec3(-2.,3.,1.), vec3(2., 1., 2.), 1);
    
    // Cylinder
    const Cylinder cy=Cylinder(vec3(0., 0., 3.), 1., 4., 1);
    
    // Ellipsoid 
    const Ellipsoid ellipsoid = Ellipsoid(vec3(-4., 0., 3.), vec3(1.5, 1.0, 0.5), 1);

    //Capsule
    const Capsule myCapsule = Capsule(vec3(3., 0.0, 3.0), vec3(5.0, 0.0, 3.0), 1., 1);

    
    x=Hit(1000.,vec3(0),-1);
    Hit current;
    bool ret=false;
    /*if(IntersectSphere(ray,sph1,current)&&current.t<x.t){
        x=current;
        ret=true;
    }*/
    
    if(IntersectSphere(ray,sph2,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if(IntersectBox(ray,bx,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if(IntersectPlane(ray,pl,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if(IntersectCylinder(ray,cy,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
    if (IntersectEllipsoid(ray, ellipsoid, current) && current.t < x.t) {
        x = current;
        ret = true;
    }
    if (IntersectCapsule(ray, myCapsule, current) && current.t < x.t) {
    x = current;
    ret = true;
}
    
    return ret;
}

vec3 Background(vec3 rd)
{
    return mix(vec3(.8,.8,.9),vec3(.7,.7,.8),rd.z);
}

// Camera rotation matrix
// ro : Camera origin
// ta : Target point
mat3 setCamera(in vec3 ro,in vec3 ta)
{
    vec3 cw=normalize(ta-ro);
    vec3 cp=vec3(0,0,1);
    vec3 cu=-normalize(cross(cw,cp));
    vec3 cv=-normalize(cross(cu,cw));
    return mat3(cu,cv,cw);
}

// Apply color model
// m : Material
// n : normal
vec3 Color(Material m,vec3 n)
{
    vec3 light=normalize(vec3(1,1,1));
    
    float diff=clamp(dot(n,light),0.,1.);
    vec3 col=m.d*diff+vec3(.2,.2,.2);
    return col;
}

// Rendering
vec3 Shade(Ray ray)
{
    // Intersect contains all the geo detection
    Hit x;
    bool idx=Intersect(ray,x);
    
    if(idx)
    {
        vec3 p=Point(ray,x.t);
        Material mat=Texture(p,x.i);
        
        return Color(mat,x.n);
    }
    else
    {
        return Background(ray.d);
    }
    
    return vec3(0);
}

void mainImage(out vec4 fragColor,in vec2 fragCoord)
{
    // From uv which are the pixel coordinates in [0,1], change to [-1,1] and apply aspect ratio
    vec2 uv=(-iResolution.xy+2.*fragCoord.xy)/iResolution.y;
    
    // Mouse control
    vec2 mouse=iMouse.xy/iResolution.xy;
    
    // Ray origin
    vec3 ro=20.*normalize(vec3(sin(2.*3.14*mouse.x),cos(2.*3.14*mouse.x),1.4*(mouse.y-.1)));
    vec3 ta=vec3(0.,0.,1.5);
    mat3 ca=setCamera(ro,ta);
    
    // Ray
    vec3 rd=ca*normalize(vec3(uv.xy*tan(radians(22.5)),1.));
    
    // Render
    vec3 col=Shade(Ray(ro,rd));
    
    fragColor=vec4(col,1.);
  

}
