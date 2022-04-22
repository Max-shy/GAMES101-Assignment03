//
// Created by LEI XU on 4/27/19.
//

#ifndef RASTERIZER_TEXTURE_H
#define RASTERIZER_TEXTURE_H
#include "global.hpp"
#include <eigen3/Eigen/Eigen>
#include <opencv2/opencv.hpp>
class Texture{
private:
    cv::Mat image_data;

public:
    Texture(const std::string& name)
    {
        image_data = cv::imread(name);
        cv::cvtColor(image_data, image_data, cv::COLOR_RGB2BGR);
        width = image_data.cols;
        height = image_data.rows;
    }

    int width, height;

    Eigen::Vector3f getColor(float u, float v)
    {
        //×ø±êÏÞ¶¨
        u = std::fmin(1, std::fmax(u, 0));
        v = std::fmin(1, std::fmax(v, 0));
        auto u_img = u * width;
        auto v_img = (1 - v) * height;
        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        return Eigen::Vector3f(color[0], color[1], color[2]);
    }

    Eigen::Vector3f Bilinear_Interpolation(float u, float v) {
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
};
#endif //RASTERIZER_TEXTURE_H
