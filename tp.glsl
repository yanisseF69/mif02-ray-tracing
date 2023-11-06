struct Sphere{
    vec3 c;// Center
    float r;// Radius
    int i;// Texture Id
};

struct Ellipsoid {
    vec3 c; // Center
    vec3 radii; // Radii along each axis (x, y, z)
    int i; // Texture Id
};

struct Plane {
    vec3 n;        // Normal
    vec3 p;        // Point
    int i;         // Texture Id
    float amplitude;
    float frequency;
    float speed;
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
    float t; // transparency
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


float random(vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898, 78.233)) * 43758.5453123));
}

float noise(vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);
    
    // Smoothstep function for smooth interpolation
    f = f * f * (3.0 - 2.0 * f);
    
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));
    
    float mix1 = mix(a, b, f.x);
    float mix2 = mix(c, d, f.x);
    
    return mix(mix1, mix2, f.y);
}

float turbulence(vec2 st, int octaves) {
    float value = 0.0;
    float amplitude = 1.0;
    for (int i = 0; i < octaves; i++) {
        value += amplitude * noise(st);
        st *= 2.0; // Increase the frequency
        amplitude *= 0.5; // Decrease the amplitude
    }
    return value;
}

float woodPattern(vec2 st) {
    // Define the scale of the wood grain
    float scale = 8.0;

    // Calculate the noise value
    float n = noise(st * scale);

    // Apply a sinusoidal function to create wood-like rings
    return abs(sin(n * 2.0 * 3.1415));
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
    } else if (i == 4) {
        // The transparent object
        res.t = 1.; 
        return res;
     } else if (i == 5) {
     // Define the scale and frequency of the marble texture
        float scale = 8.0;
        int octaves = 5;

        // Calculate the turbulence value
        float t = turbulence(p.xy * scale, octaves);

        // Define the marble colors
        vec3 darkColor = vec3(0.1, 0.1, 0.1);
        vec3 lightColor = vec3(1.0, 1.0, 1.0);

        // Mix the colors based on the turbulence
        res.d = mix(darkColor, lightColor, t);
        return res;
     
     } else if (i == 6) {
     // Generate a wood grain pattern
        float woodGrain = woodPattern(p.xy);

        // Define the wood colors
        vec3 lightColor = vec3(0.6, 0.3, 0.1);
        vec3 darkColor = vec3(0.3, 0.1, 0.05);

        // Mix the colors based on the wood grain pattern
        res.d = mix(darkColor, lightColor, woodGrain);
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
bool IntersectPlane(Ray ray, Plane pl, out Hit x) {
    float t = -dot(ray.o - pl.p, pl.n) / dot(ray.d, pl.n);
    if (t > 0. && t < 1000.0) {
        // Introduce waves that move with time
        float waveHeight = pl.amplitude * sin(pl.frequency * t + pl.speed * iTime);

        // Modify the plane's position based on the wave
        vec3 wavePosition = pl.p + pl.n * waveHeight;

        x = Hit(t, pl.n, pl.i);
        x.n = pl.n;
        x.n.y = waveHeight;  // Adjust the normal to account for the wave
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
    vec3 oc = ray.o - ellipsoid.c;  // Translate ray origin
    oc /= ellipsoid.radii;  // Scale by radii

    vec3 d = ray.d / ellipsoid.radii;  // Scale ray direction

    float a = dot(d, d);
    float b = 2.0 * dot(oc, d);
    float c = dot(oc, oc) - 1.0;

    float discriminant = b * b - 4.0 * a * c;

    if (discriminant > 0.0) {
        float t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        float t2 = (-b + sqrt(discriminant)) / (2.0 * a);

        if (t1 > 0.0 || t2 > 0.0) {
            float t = (t1 > 0.0) ? t1 : t2;
            vec3 p = ray.o + t * ray.d;
            vec3 normal = normalize(2.0 * (p - ellipsoid.c));

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

        //if (intersectZLocal > -cy.h && intersectZLocal < cy.h ) {
            vec3 p = Point(ray, t);
            h = Hit(t, normalize(p - cy.c), cy.i);
            return true;
        //}
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
    Plane pl = Plane(vec3(0., 0., 1.), vec3(0., 0., 0.), 3, 0.25, 5.0, 2.5);

    
    // Sphere
    Sphere s = Sphere(vec3(7., -7., 7), 3., 1);

    // Box
    Box bx = Box(vec3(3., 6., 2.), vec3(7., 3., 4.), 1);

    // Cylinder
    const Cylinder cy = Cylinder(vec3(0., 0., 3.), 1., 4., 6);

    // Ellipsoid 
    const Ellipsoid ellipsoid = Ellipsoid(vec3(-5., 0., 7.5), vec3(3., 2., 2.), 1);

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
    vec2 translation = currentPoint - points[0]; 
    bx.mini.xy += translation;
    bx.maxi.xy += translation;

    
    // Torus homothétie
    float rotationAngleCapsule = iTime;
    mat2 rotationMatrixCapsule = mat2(
        cos(rotationAngleCapsule), -sin(rotationAngleCapsule),
        sin(rotationAngleCapsule), cos(rotationAngleCapsule)
    );

    vec2 rotatedA = rotationMatrixCapsule * capsule.a.xy;
    vec2 rotatedB = rotationMatrixCapsule * capsule.b.xy;

    capsule.a.xyz = vec3(rotatedA.x, rotatedA.y, capsule.a.z);
    capsule.b.xyz = vec3(rotatedB.x, rotatedB.y, capsule.b.z);



    if (IntersectSphere(ray, s, current) && current.t < x.t) {
        x = current;
        ret = true;
    }
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


// Hemisphere direction
vec3 Hemisphere(int seed,vec3 n)

{

    float a=fract(sin(176.19*float(seed)));// Uniform randoms
    float b=fract(sin(164.19*float(seed)));

    float u=2.*3.1415*a;// Random angle
    float v=acos(2.*b-1.);// Arcosine distribution to compensate for poles

    vec3 d=vec3(cos(u)*cos(v),sin(u)*cos(v),sin(v));// Direction

    if(dot(d,n)<0.){d=-d;}// Hemishpere

    return d;

}

// Ambient occlusion
// p : Point
// n : Normal
// N : Number of samples
float AmbientOcclusion(vec3 p, vec3 n, int N)
{
    if (N == 0) {
        return 1.0;
    }

    float ao = 0.0;

    for (int i = 0; i < N; i++)
    {
        // Generate a random direction in the hemisphere
        vec3 d = Hemisphere(i, n);

        // Create a shadow ray from the point of intersection with an offset along the normal
        vec3 offsetPoint = p + n;
        Ray shadowRay;
        shadowRay.o = offsetPoint;
        shadowRay.d = d;

        // Check if the shadow ray intersects any objects
        Hit shadowHit;
        bool shadow = Intersect(shadowRay, shadowHit);

        // If there is no intersection or the intersection is far enough, accumulate occlusion
        if (!shadow || shadowHit.t > length(p - shadowRay.o)) {
            ao += 1.0;
        }
    }

    // Normalize the accumulated occlusion factor
    ao /= float(N);

    return ao;
}





// Rendering
vec3 Shade(Ray ray) {
    int maxDepth = 2; // Maximum reflection depth
    vec3 color = vec3(0.0);

    for (int depth = 0; depth < maxDepth; depth++) {
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

            bool shadow = (Intersect(shadowRay, shadowHit) && shadowHit.t > length(lightDir));

            if (shadow) {
                // In shadow, return ambient color only
                color += mat.a;
                break; // Exit the loop if in shadow
            }

            if (mat.t > 0.0) {
                // Calculate the reflected ray direction
                vec3 reflectionDir = reflect(ray.d, x.n);

                // Create a reflected ray
                Ray reflectedRay;
                reflectedRay.o = p + 0.001 * reflectionDir; // Offset to avoid self-intersection
                reflectedRay.d = reflectionDir;

                // Trace the reflected ray and accumulate the color
                color += mat.d * (1.0 - mat.t);
                ray = reflectedRay; // Update the ray for the next iteration
            } else {
                // Calculate object's color with ambient occlusion
                float ao = AmbientOcclusion(p, x.n, 0); // Adjust the number of samples as needed
                // Apply the occlusion factor to the color
                color += Color(mat, x.n, p) * ao;
                break; // Exit the loop if not transparent
            }
        } else {
            // No intersection, return the background color
            color += Background(ray.d);
            break; // Exit the loop if no intersection
        }
    }

    return color;
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

