# GAMES101-Homework03
实现Blinn-Phong模型进行着色，并进行光栅化

Assignment 3 requires the interpolation of vectors, colors, and textures during rasterization.  I need to complete the Blinn-Phong reflection model and modify shaders to achieve different texture effects. 

I don't need to mention the perspective transformation part, just copy the code here. Now, I need to modify the function **rasterize_triangle()** to interpolate colors, normal vectors, textures, and texture coordinates. 

```CPP
//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t, const std::array<Eigen::Vector3f, 3>& view_pos) 
{
    // TODO: From your HW3, get the triangle rasterization code.
    // TODO: Inside your rasterization loop:
    auto v = t.toVector4();
    //create the Bounding Box
    float x_min = std::min(t.v[0].x(), std::min(t.v[1].x(), t.v[2].x()));
    float x_max = std::max(t.v[0].x(), std::max(t.v[1].x(), t.v[2].x()));
    float y_min = std::min(t.v[0].y(), std::min(t.v[1].y(), t.v[2].y()));
    float y_max = std::max(t.v[0].y(), std::max(t.v[1].y(), t.v[2].y()));
    int xmin = std::floor(x_min);
    int xmax = std::ceil(x_max);
    int ymin = std::floor(y_min);
    int ymax = std::ceil(y_max);

    for (int i = xmin; i <= xmax; i++) {
        for (int j = ymin; j <= ymax; j++) {
            if (insideTriangle(i + 0.5, j + 0.5, t.v)) {
                 auto[alpha, beta, gamma] =computeBarycentric2D(i+0.5,j+0.5,t.v);
                 float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                 float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                 z_interpolated *= w_reciprocal; //z-buffer

                 if (z_interpolated < depth_buf[get_index(i, j)])
                 {      
                     //interpolate
                     // color 
                     auto interpolated_color = interpolate(alpha, beta, gamma, t.color[0], t.color[1], t.color[2], 1);
                     // normal
                     auto interpolated_normal = interpolate(alpha, beta, gamma, t.normal[0], t.normal[1], t.normal[2], 1).normalized();
                     // texture
                     auto interpolated_texcoords = interpolate(alpha, beta, gamma, t.tex_coords[0], t.tex_coords[1], t.tex_coords[2], 1);
                     // texture coords
                     auto interpolated_shadingcoords = interpolate(alpha, beta, gamma, view_pos[0], view_pos[1], view_pos[2], 1);

                     // Use a structure to Save interpolation data,
                     fragment_shader_payload payload(interpolated_color, interpolated_normal, interpolated_texcoords, texture ? &*texture : nullptr);

                     payload.view_pos = interpolated_shadingcoords;
                     auto pixel_value = fragment_shader(payload); //Store the various attributes of the pixel in this variable
                     // set the z-buffer
                     depth_buf[get_index(i, j)] = z_interpolated;
                     // set the pixel various,(colors, normal, textures, texture coordinates)
                     set_pixel(Eigen::Vector3f(i, j, 0), pixel_value);
                 }
            }
        }
    }
}
```

I didn't use MSAA here because it was too slow and not very effective.

**normal_shader_shader()**:

![output_normal_fragment_shader](https://user-images.githubusercontent.com/68177870/168940869-2d97fef1-15c2-45ff-8396-f20201f6dc40.png)



In the next part, I need to complete the Blinn-Phong reflection model in the function **phong_fragment_shader()**.

I'll start with a review of the Blinn-Phong reflection model.

**Blinn-Phong reflection**

It consists of three parts:

- Diffuse reflection: To reflect light around with the same intensity
 
- Specular highlights: The brightest part
  
- Ambient light: A constant color of light

The result color is equal to the sum of these three parts.


```CPP
Eigen::Vector3f phong_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);//Ambient light parameter
    Eigen::Vector3f kd = payload.color / 255.f;//diffuse reflection parameter
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);//specular light parameter

    auto l1 = light{{20, 20, 20}, {500, 500, 500}};
    auto l2 = light{{-20, 20, 0}, {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};// Ia
    Eigen::Vector3f eye_pos{0, 0, 10};//camera position
    float p = 150;
    Eigen::Vector3f color = payload.color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;
    Eigen::Vector3f result_color = {0, 0, 0};
    for (auto& light : lights)
    {
        // TODO: For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular* 
        // components are. Then, accumulate that result on the *result_color* object.
        Eigen::Vector3f light_dir = light.position - point;//l,light direction
        Eigen::Vector3f view_dir = eye_pos - point;//v , camera direction
        float r = light_dir.dot(light_dir);// r^2,distance of light and point 

        //ambient light
        Eigen::Vector3f La = ka.cwiseProduct(amb_light_intensity);//Dot product
        //diffuse
        Eigen::Vector3f Ld = kd.cwiseProduct(light.intensity / r);
        Ld *= std::max(0.0f, normal.normalized().dot(light_dir.normalized()));
        //speculor 
        Eigen::Vector3f h = (light_dir + view_dir).normalized();//add two vector addition, get angle bisector
        Eigen::Vector3f Ls = ks.cwiseProduct(light.intensity / r);
        Ls *= std::pow(std::max(0.0f, normal.normalized().dot(h)), p);
        result_color += (La + Ld + Ls);
    }
    return result_color * 255.f;
}
```

phong_fragment_shader result: 

![output_phong_fragment_shader](https://user-images.githubusercontent.com/68177870/168940969-57c526ae-319a-4e1f-82f8-496254681c6b.png)





Next, I need to load the texture and let parameter kd = texture_color in the function **texture_fragment_shader().**

```CPP
Eigen::Vector3f texture_fragment_shader(const fragment_shader_payload& payload)
{
    Eigen::Vector3f return_color = { 0, 0, 0 };
    if (payload.texture)
    {
        // TODO: Get the texture value at the texture coordinates of the current fragment
        return_color = payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y());
    }
    Eigen::Vector3f texture_color;
    texture_color << return_color.x(), return_color.y(), return_color.z();

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);//Ambient light parameter
    Eigen::Vector3f kd = texture_color / 255.f;//diffuse reflection parameter
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);//specular light parameter

   //....same as phong_fragment_shader

    for (auto& light : lights)
    {
        //....same asphong_fragment_shader
        //calculate the *ambient*, *diffuse*, and *specular* 
    }
    return result_color * 255.0f;
}

```

texture_fragment_shader result: 

![output_texture_fragment_shader](https://user-images.githubusercontent.com/68177870/168940970-095bbf6b-4f97-4cf5-a710-0f2d4a3a8afa.png)





Next, I need to complete the rendering with the **Bump mapping**.

Bump mapping makes the object concave and convex by adding disturbance to the normal vector.

Set the initial normal vector to  (0,0,1), according to the calculation method in 3D: 

The updated normal vector is : 

Then use the  inverse  view matrix transforms the normal vector to the global coordinates: 

And then just normalize the normal vectors.

```CPP
Eigen::Vector3f bump_fragment_shader(const fragment_shader_payload& payload)
{
    //....same as phong_fragment_shader
	//set the l1,l2,ka,ks,kd.....
    float kh = 0.2, kn = 0.1;

    // TODO: Implement bump mapping here
    //adding disturbance to the normal vector
    float u = payload.tex_coords.x();
    float v = payload.tex_coords.y();
    float w = payload.texture->width;
    float h = payload.texture->height;

    float du = kh * kn * (payload.texture->getColor(u + 1.0f / w, v).norm() - payload.texture->getColor(u, v).norm());
    float dv = kh * kn * (payload.texture->getColor(u, v + 1.0f / h).norm() - payload.texture->getColor(u, v).norm());
    Eigen::Vector3f ln(-du, dv, 1); //the new normal

    float x = normal.x(),y = normal.y(),z = normal.z();
    Eigen::Vector3f g_t(x * y / std::sqrt(x * x + z * z),std::sqrt(x * x + z * z),z * y / std::sqrt(x * x + z * z));
    Eigen::Vector3f t = normal.cross(g_t);

    //view Matrix，set the vector (0,0,1) to world coordinates
    Eigen::Matrix3f TEG;
    TEG << g_t.x(), t.x(), normal.x(),
           g_t.y(), t.y(), normal.y(),
           g_t.z(), t.z(), normal.z();
    normal = TEG * ln;

    result_color = normal.normalized();

    return result_color * 255.f;
}
```

bump_fragment_shader result:

![output_bump_fragment_shader](https://user-images.githubusercontent.com/68177870/168940990-b052b24f-e116-4488-ad77-5c56ca3ea482.png)




Displacement mapping is a more advanced approach.

- Actually moves the vertices.
- Use the same texture as in bumping mapping.

```CPP
Eigen::Vector3f displacement_fragment_shader(const fragment_shader_payload& payload)
{
    
    //....same as bump_fragment_shader
	//adding disturbance to the normal vector
    
    //Actually moves the point vertices.
    point += (kn * normal* payload.texture->getColor(u, v).norm());
    //view Matrix，set the vector (0,0,1) to world coordinates.
    normal = (TEG * ln).normalized();

    Eigen::Vector3f result_color = { 0, 0, 0 };
    for (auto& light : lights)
    {
        //....same asphong_fragment_shader
        //calculate the *ambient*, *diffuse*, and *specular* 
    }
    return result_color * 255.f;
}
```

Displacement mapping result:

![output_displacement_fragment_shader](https://user-images.githubusercontent.com/68177870/168940995-2656002c-b1f4-489e-bc4b-ecf5fe6421ee.png)




I tried to render the other models offered in the course. 

1. rendering bunny.obj with normal_fragment_shader.

  ![output](https://user-images.githubusercontent.com/68177870/168941029-38afc78b-434d-44dc-a175-810453a2d002.png)


2. rendering spot control_mesh.obj with normal_fragment_shader.

   

   

   There were problems when I rendered other models(rock.obj). It seems that the model vertices do not correspond to texture coordinates.



At last, I tried to implement Bilinear Interpolation in the Texture class.


We need to take four pixels and interpolate the pixel values within them.

```CPP
 Eigen::Vector3f BilinearInterpolation(float u, float v) {
        //takes four pixel 
        float u01 =int(u*width),v01 = int(v*height) ;
        float u11 = u01 + 1, v11 = v01;
        float u00 = u01, v00 = v01 + 1;
        float u10 = u01 + 1, v10 = v01 + 1;

        //Calculate the color of four pixels
        Eigen::Vector3f color01 = getColor(u01 / width, v01 / height);
        Eigen::Vector3f color11 = getColor(u11 / width, v11 / height);
        Eigen::Vector3f color00 = getColor(u00 / width, v00 / height);
        Eigen::Vector3f color10 = getColor(u10 / width, v10 / height);

        //Bilinear Interpolation
        Eigen::Vector3f Coloru1 = color01 + (color11 - color01) * (u * width - u01);
        Eigen::Vector3f Coloru0 = color00 + (color10 - color00) * (u * width - u00);
        Eigen::Vector3f Color_result = Coloru1 + (Coloru0 - Coloru1) * (v * height - v01);
        return Color_result;
    }
```

Compare bilinear Interpolation effects with different textures :

before:





After:




We can definitely feel the image becoming clearer.





Using bump_fragment_shader to compare bilinear Interpolation effects :

before:

![output_bump_fragment_shader](E:\CG\Games\GAMES101\Homework\03(Lecture 06)\GAMES101-Homework03\images\spot\output_bump_fragment_shader.png)



After:

![output_Bilinear_Interpolation2](E:\CG\Games\GAMES101\Homework\03(Lecture 06)\GAMES101-Homework03\images\spot\output_Bilinear_Interpolation2.png)

We can definitely feel the image becoming smoother.



