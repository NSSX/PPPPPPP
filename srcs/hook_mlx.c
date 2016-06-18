/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   hook_mlx.c                                         :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/05 14:17:05 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 14:17:16 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include "mlx_encoding.h"
#include "mlx.h"
#include <stdio.h>
#include <math.h>

static void	ft_apply_gamma(t_image *image, int x, int y)
{
	return ;
	int			offset;
	t_ucolor	uc;
	t_fcolor	fc;

	if (x < 0 || y < 0 || x >= image->width || y >= image->height)
		return ;
	offset = (y * image->width + x) * image->bypp;
	uc.b = (image->data + offset)[0];
	uc.g = (image->data + offset)[1];
	uc.r = (image->data + offset)[2];
	fc.r = (float)uc.r / 255.0;
	fc.b = (float)uc.b / 255.0;
	fc.g = (float)uc.g / 255.0;
	(image->data + offset)[0] = powf(fc.b, 1.9) * 255.0;
	(image->data + offset)[1] = powf(fc.g, 1.9) * 255.0;
	(image->data + offset)[2] = powf(fc.r, 1.9) * 255.0;
}

/*

	int row = blockIdx.y * blockDim.y + threadIdx.y;
  	int col = blockIdx.x * blockDim.x + threadIdx.x;

  	if(col >= WIDTH || row >= HEIGHT) return;

  	int imindex;
  	imindex = row * linesize;
	imindex += (col * bpx / 8);

  	float rweight;
	float gweight;
	float bweight;
	unsigned char r = addr[imindex];
	unsigned char g = addr[imindex + 1];
	unsigned char b = addr[imindex + 2];

	int gray = 0;
	if (gray_scale_method == AVERAGE)
	{

		gray = (r + g + b)/3;
	}
	else if (gray_scale_method == GIMP)
	{
		rweight = 0.3;
		gweight = 0.59;
		bweight = 0.11;
		gray = (r * rweight + g * gweight + b * bweight);
	}
	else if (gray_scale_method == LUMA)
	{
		rweight =0.2126;
		gweight = 0.7152;
		bweight = 0.0722;
		gray = (r * rweight + g * gweight + b * bweight);
	}
	else if (gray_scale_method == BT601)
	{
		rweight =0.299;
		gweight = 0.587;
		bweight = 0.114;
		gray = (r * rweight + g * gweight + b * bweight);
	}
	else if (gray_scale_method == SHADE)
	{
		float conversion_factor = 255 / (nshade - 1);
		float average = (r + g + b) / 3;
		gray = (int)((average / conversion_factor) + 0.5) * conversion_factor;
	}
	addr[imindex] = gray;
	addr[++imindex] = gray;
	addr[++imindex] = gray;



*/

static void	ft_post_processing(t_env *e)
{
	int x;
	int y;

	for (x = 0; x < e->screen->width; ++x)
	{
		for (y = 0; y < e->screen->height; ++y)
		{
			ft_apply_gamma(e->screen, x, y);
		}
	}
}

int		expose_hook(t_env *e)
{
	mlx_put_image_to_window(e->mlx, e->win, e->screen->img, 0, 0);
	return (0);
}

int		loop_hook(t_env *e)
{
	if (FLAG_GET(e->flag, MASK_RECALC))
	{
		recalc_scene(e);
		ft_post_processing(e);
		ft_printelapsed(e);
		FLAG_UNSET(e->flag, MASK_RECALC);
	}
	mlx_put_image_to_window(e->mlx, e->win, e->screen->img, 0, 0);
	return (0);
}

int		key_hook(int keycode, t_env *e)
{
	if (keycode == KEYCODE_ESCAPE)
		ft_exit(0, e);
	return (0);
}
