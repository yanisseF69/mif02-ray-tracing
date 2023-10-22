struct Ellipsoid {
    vec3 c; // Center
    vec3 radii; // Radii along each axis (x, y, z)
    int i; // Texture Id
};

struct Plane {
    vec3 n; // Normal
    vec3 p; // Point
    int i; // Texture Id
};

struct Hit {
    float t; // Intersection depth
    vec3 n; // Normal
    int i; // Texture Id
};

struct Ray {
    vec3 o; // Origin
    vec3 d; // Direction
};

struct Material {
    vec3 a, d, s; 
};

struct Box {
    vec3 mini;
    vec3 maxi;
    int i;
};

struct Cylinder {
    vec3 c;
    float r;
    float h;
    int i;
};

struct Torus {
    vec3 c;
    vec2 t;
    int i;
};

struct Capsule {
    vec3 a; // debut 
    vec3 b; // fin 
    float r; // radius
    int i; // texture id
};

float Checkers(in vec2 p) {
    // Filter kernel
    vec2 w = fwidth(p) + .001;
    // Box box filter
    vec2 i = 2. * (abs(fract((p - .5 * w) * .5) - .5) - abs(fract((p + .5 * w) * .5) - .5)) / w;
    // xor pattern
    return .5 - .5 * i.x * i.y;
}

// Compute point on ray
vec3 Point(Ray ray, float t) {
    return ray.o + t * ray.d;
}

int convert(float x)
{
    if(x < 0.0) x--;
    return int(x);
}

// Compute color
// i : Texture index
// p : Point
Material Texture(Ray ray, vec3 p, int i) {
    
    Material res;
    res.a = vec3(.1, .1, .1);
    res.s = vec3(1.0, 1, 1.);

    if (i == 1) {
        res.d = vec3(.8, .5, .4);
        return res;
    } else if (i == 0) {

        float f = Checkers(.5 * p.xy);
        Material mat;
        res.d = vec3(.4, .5, .7) + f * vec3(.1);
        return res;
    } else if (i == 2) {
        
        vec3 noir = vec3(0., 0., 0.);
        vec3 blanc = vec3(1., 1., 1.);

        int x = convert(p.x);
        int y = convert(p.y);
        int z = convert(p.z);

        if ((x + y + z) % 2 == 0) {
            res.d = noir;
            return res;
        } else {
        
            res.d = blanc; // Couleur diffuse
            return res;
            }
    } else if (i == 3) {
        // Compute water effect with deformations
        vec2 uv = 0.5 * p.xy; // Scale down the input coordinates

        // Deform the UV coordinates with sine waves
        float time = iTime; // Assuming you have a uniform called iTime for animation
        uv.x += 0.1 * sin(uv.y * 10.0 + time * 2.0);
        uv.y += 0.1 * sin(uv.x * 10.0 + time * 2.0);

        res.d = vec3(0.0, 0.3, 0.7) + 0.1 * sin(uv.x * 10.0) + 0.1 * sin(uv.y * 10.0);
        return res;
    } else {
        // Sample from the uniform sampler iChannel0 using the texture coordinates
        vec2 texCoord = p.xy;
        float f = texture(iChannel0, texCoord).r;

        res.d = vec3(.4, .5, .7) + f * vec3(.1);
        return res;
    }
    
    res.d = vec3(0.0);
    return res;
}

// Plane intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectPlane(Ray ray, Plane pl, out Hit x) {
    float t = -dot(ray.o - pl.p, pl.n) / dot(ray.d, pl.n);
    if (t > 0.) {
        x = Hit(t, vec3(0, 0, 1), 0);
        return true;
    }
    return false;
}

bool IntersectBox(Ray ray, Box box, out Hit x) {

    float tMin = -99999.0; //distance d entree la plus proche
    float tMax = 99999.0; //distance de sortie la plus lointaine

    //on calcule les potentielles distances d'entrées et de sortie
    //pour chaque dimensions
    for (int i = 0; i < 3; i++) {

        float t0 = (box.mini[i] - ray.o[i]) / ray.d[i];
        float t1 = (box.maxi[i] - ray.o[i]) / ray.d[i];

        if (t1 < t0) {

            float tmp = t1;
            t1 = t0;
            t0 = tmp;
        }

        tMin = max(tMin, t0);
        tMax = min(tMax, t1);

        if (tMin > tMax) {

            return false;
        }
    }

    vec3 p = Point(ray, tMin);
    x = Hit(tMin, normalize(p - (box.mini + box.maxi) * 0.5), box.i);

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

bool IntersectCylinder(Ray ray, Cylinder cy, out Hit h) {
    vec3 rayOriginLocal = ray.o - cy.c;

    float a = ray.d.x * ray.d.x + ray.d.y * ray.d.y;
    float b = 2.0 * (ray.d.x * rayOriginLocal.x + ray.d.y * rayOriginLocal.y);
    float c = rayOriginLocal.x * rayOriginLocal.x + rayOriginLocal.y * rayOriginLocal.y - cy.r * cy.r;

    float delta = b * b - 4.0 * a * c;

    if (delta > 0.0) {
        float t1 = (-b - sqrt(delta)) / (2.0 * a);
        float t2 = (-b + sqrt(delta)) / (2.0 * a);
        float t;

        if (t1 < t2) t = t1;
        else t = t2;

        float intersectXLocal = ray.d.x * t + rayOriginLocal.x;
        float intersectYLocal = ray.d.y * t + rayOriginLocal.y;
        float intersectZLocal = ray.d.z * t + rayOriginLocal.z;

        if (intersectXLocal >= -cy.h / 2.0 && intersectXLocal <= cy.h / 2.0) {
            vec3 p = Point(ray, t);
            h = Hit(t, normalize(p - cy.c), cy.i);
            return true;
        }
    }

    return false;
}

bool IntersectTorus(Ray ray, Torus tor, out Hit x) {
    float tMin = 0.0;
    float tMax = 1000.0;  // une distance max arbitraire pour l'intersection
    float epsilon = 0.01; // une petite valeur pour vérifier la proximité à la surface

    for (int i = 0; i < 100; ++i) {  // Iterations limitées pour éviter une boucle infinie
        vec3 p = Point(ray, tMin);
        vec3 localP = p - tor.c;

        // This is essentially the body of the sdTorus function inlined
        vec2 q = vec2(length(localP.xz) - tor.t.x, localP.y);
        float dist = length(q) - tor.t.y;

        if (dist < epsilon) {
            // Calculate normal using finite differences
            vec3 nor = normalize(vec3(
                length(vec2(length(vec3(localP.x + epsilon, localP.y, localP.z).xz) - tor.t.x, vec3(localP.x + epsilon, localP.y, localP.z).y)) - tor.t.y - dist,
                length(vec2(length(vec3(localP.x, localP.y + epsilon, localP.z).xz) - tor.t.x, vec3(localP.x, localP.y + epsilon, localP.z).y)) - tor.t.y - dist,
                length(vec2(length(vec3(localP.x, localP.y, localP.z + epsilon).xz) - tor.t.x, vec3(localP.x, localP.y, localP.z + epsilon).y)) - tor.t.y - dist
            ));

            x = Hit(tMin, nor, tor.i);
            return true;
        }

        tMin += dist;

        if (tMin >= tMax) {
            return false;
        }
    }

    return false;
}

// Capsule intersection
bool IntersectCapsule(Ray ray, Capsule cap, out Hit x) {
    float tMin = 0.0;
    float tMax = 1000.0; // une distance max arbitraire pour l'intersection
    float epsilon = 0.01; // une petite valeur pour vérifier la proximité à la surface

    for(int i = 0; i < 100; ++i) { // Iterations limitées pour éviter une boucle infinie
        vec3 p = Point(ray, tMin);
        vec3 pa = p - cap.a;
        vec3 ba = cap.b - cap.a;
        float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
        float dist = length(pa - ba * h) - cap.r;

        if (dist < epsilon) {
            x = Hit(tMin, normalize(p - (cap.a + ba * h)), cap.i);
            return true;
        }

        tMin += dist;

        if (tMin >= tMax) {
            return false;
        }
    }
    
    return false;
}

// Scene intersection
// ray : The ray
//   x : Returned intersection information
bool Intersect(Ray ray, out Hit x) {


    // Plane
    const Plane pl = Plane(vec3(0., 0., 1.), vec3(0., 0., 0.), 0);

    // Box
    Box bx = Box(vec3(3., 6., 1.), vec3(7., 3., 2.), 1);

    // Cylinder
    const Cylinder cy = Cylinder(vec3(0., 0., 3.), 1., 4., 2);

    // Ellipsoid 
    const Ellipsoid ellipsoid = Ellipsoid(vec3(-3., 0., 3.5), vec3(1.5, 1.0, 0.5), 1);

    // Torus 
    Torus torus = Torus(vec3(3., 0., 7.0), abs(vec2(1.0, 0.5)*cos(iTime)), 1);

    // Capsule 
    Capsule capsule = Capsule(vec3(-3., 0., 5.), vec3(-6., 0., 5.), 0.5, 1);

    x = Hit(1000., vec3(0), -1);
    bool ret = false;
    Hit current;


    if (IntersectPlane(ray, pl, current) && current.t < x.t) {
        x = current;
        ret = true;
    }

    // Define an array of 4 points
    vec2 points[4];
    points[0] = vec2(5.0, 5.0);
    points[1] = vec2(5.0, -5.0);
    points[2] = vec2(-5.0, -5.0);
    points[3] = vec2(-5.0, 5.0);

    // Calculate the current index and interpolation factor based on time
    float cycleTime = -4.0; // Time to complete one cycle
    float t = mod(iTime, cycleTime) / cycleTime;
    int currentIndex = int(floor(t * 4.0));
    float interpFactor = fract(t * 4.0);

    // Interpolate between the current and next points
    vec2 currentPoint = mix(points[currentIndex], points[(currentIndex + 1) % 4], interpFactor);

    // Apply translation to the box vertices
    vec2 translation = currentPoint - points[0]; // Calculate the translation relative to the first point
    bx.mini.xy += translation;
    bx.maxi.xy += translation;
    ////////////////////////////////////////////////////////////////////
    
    // Calculate rotation angle for sph1 around the Z-axis
    /*float rotationAngleSphere = iTime; // Adjust the rotation speed as needed
    mat3 rotationMatrixSphere = mat3(
        cos(rotationAngleSphere), -sin(rotationAngleSphere), 0.0,
        sin(rotationAngleSphere), cos(rotationAngleSphere), 0.0,
        0.0, 0.0, 1.0
    );

    sph2.c = rotationMatrixSphere * sph2.c; // Apply rotation to sph1 position */
    
    ////////////////////////////////////////////////////////////////////
    
    // Torus homothétie

    
    // Calculate rotation angle for the capsule around the Z-axis
    float rotationAngleCapsule = iTime; // Adjust the rotation speed as needed
    mat2 rotationMatrixCapsule = mat2(
        cos(rotationAngleCapsule), -sin(rotationAngleCapsule),
        sin(rotationAngleCapsule), cos(rotationAngleCapsule)
    );

    vec2 rotatedA = rotationMatrixCapsule * capsule.a.xy;
    vec2 rotatedB = rotationMatrixCapsule * capsule.b.xy;

    capsule.a.xyz = vec3(rotatedA.x, rotatedA.y, capsule.a.z);
    capsule.b.xyz = vec3(rotatedB.x, rotatedB.y, capsule.b.z);




    if (IntersectBox(ray, bx, current) && current.t < x.t) {
        x = current;
        ret = true;
    }

    if (IntersectCylinder(ray, cy, current) && current.t < x.t) {
        x = current;
        ret = true;
    }

    if (IntersectEllipsoid(ray, ellipsoid, current) && current.t < x.t) {
        x = current;
        ret = true;
    }

    if (IntersectTorus(ray, torus, current) && current.t < x.t) {
        x = current;
        ret = true;
    }
    if (IntersectCapsule(ray, capsule, current) && current.t < x.t) {
        x = current;
        ret = true;
    }

    return ret;
}

vec3 Background(vec3 rd) {
    return mix(vec3(.8, .8, .9), vec3(.7, .7, .8), rd.z);
}

// Camera rotation matrix
// ro : Camera origin
// ta : Target point
mat3 setCamera(in vec3 ro, in vec3 ta) {
    vec3 cw = normalize(ta - ro);
    vec3 cp = vec3(0, 0, 1);
    vec3 cu = -normalize(cross(cw, cp));
    vec3 cv = -normalize(cross(cu, cw));
    return mat3(cu, cv, cw);
}

// Apply color model
// m : Material
// n : normal
vec3 Color(Material m, vec3 n, vec3 v)
{
    vec3 lightDir = normalize(vec3(1, 1, 1));
    vec3 viewDir = normalize(v);
    vec3 reflectionDir = reflect(-lightDir, n);
    
    float diff = max(dot(n, lightDir), 0.0);
    float spec = 0.0;
    
    if (diff > 0.0) {
        spec = pow(max(dot(viewDir, reflectionDir), 0.0), 32.0);
    }
    
    vec3 ambient = m.a;
    vec3 diffuse = m.d * diff;
    vec3 specular = m.s * spec;

    vec3 col = ambient + diffuse + specular;
    
    return col;
}

// Rendering
vec3 Shade(Ray ray) {
    // Intersect contains all the geo detection
    Hit x;
    bool idx = Intersect(ray, x);

    if (idx) {
        vec3 p = Point(ray, x.t);
        Material mat = Texture(ray, p, x.i);
        vec3 lightDir = normalize(vec3(1., 1., 1.));
        Ray shadowRay;
        shadowRay.o = p + 0.001 * lightDir; // Offset the origin slightly to avoid self-intersection
        shadowRay.d = lightDir;

        Hit shadowHit;
        
        bool shadow = (Intersect(shadowRay, shadowHit) && shadowHit.t >= length(lightDir));
        
        if (shadow) {
            // In shadow, return ambient color only
            return mat.a;
        } else {
            // Not in shadow, calculate and return object's color
            return Color(mat, x.n, p);
        }
    } else {
        return Background(ray.d);
    }

    return vec3(0);
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // From uv which are the pixel coordinates in [0,1], change to [-1,1] and apply aspect ratio
    vec2 uv = (-iResolution.xy + 2. * fragCoord.xy) / iResolution.y;

    // Mouse control
    vec2 mouse = iMouse.xy / iResolution.xy;

    // Ray origin
    vec3 ro = 25. * normalize(vec3(sin(2. * 3.14 * mouse.x), cos(2. * 3.14 * mouse.x), 1.4 * (mouse.y - .1)));
    vec3 ta = vec3(0., 0., 1.5);
    mat3 ca = setCamera(ro, ta);

    // Ray
    vec3 rd = ca * normalize(vec3(uv.xy * tan(radians(22.5)), 1.));

    // Render
    vec3 col = Shade(Ray(ro, rd));

    fragColor = vec4(col, 1.);
}
