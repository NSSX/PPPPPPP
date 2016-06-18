/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   rtv1.h                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/05 14:31:48 by igomez            #+#    #+#             */
/*   Updated: 2016/06/17 15:18:59 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef RTV1_H
# define RTV1_H

# define TITLE 					"RTv1"
# define WIN_WIDTH				500
# define WIN_HEIGHT 			WIN_WIDTH
# define DIST_TO_PROJ			WIN_WIDTH / 2

# define FLAG_GET(F, M)			(!!(F & M))
# define FLAG_SET(F, M)			(F |= M)
# define FLAG_UNSET(F, M)		(F &= ~M)
# define FLAG_TOGGLE(F, M)		(F ^= M)
# define FLAG_PRINT(F)			ft_flag_print(&F)
# define MASK_RECALC			(1U << 0)

# define ERROR					0.000000001
# define COLOR_BACKGROUND		0x0
# define INTENSITY_MIN			0.1

# include <time.h>
# include <stdint.h>
# include "libft.h"
# include "rtv1_matrix.h"
# include "rtv1_struct.h"

/*
** Env structure
*/

typedef struct		s_env
{
	double			nb_reflexion;
	int				width;
	int				height;
	char			*title;
	void			*mlx;
	void			*win;
	t_image			*screen;
	int				flag;
	clock_t			c_old;
	clock_t			c_new;
	t_scene			*scene;
}					t_env;

/*
** Initialisation
*/

t_image				*image_new(t_env *e, int width, int height);
int					env_mlx_init(t_env *e);
int					env_scene_init(t_env *e);

/*
** Cleaning and exiting functions
*/

void				clean_env(t_env *e);
void				ft_exit(int i, t_env *e);

/*
** MLX Events functions
*/

int					expose_hook(t_env *e);
int					key_hook(int keycode, t_env *e);
int					loop_hook(t_env *e);

/*
** Drawing functions
*/

int					recalc_scene(t_env	*e);

/*
** Ray functions
*/

int					ray_setup_camera(t_ray *ray, t_cam *cam);
int					ray_setup_obj(t_ray *ray, t_obj *obj);
int					ray_cleanup_obj(t_ray *ray, t_obj *obj);
int					ray_intersect(t_ray *ray, t_env *e);
int					ray_color(t_ray *ray, t_env *e);
double				solve_2_degree(double a, double b, double c);

/*
** Color
*/

int					color_scale(int c, double intensity);
int					color_add(int a, int b);

/*
** Debug only
*/

void				ft_printmat(t_matrix *buf);
void				ft_printobj(t_obj *obj);
void				ft_printray(t_ray *ray);
void				ft_printelapsed(t_env *e);
void				flat_color(int x, int y, t_obj *obj, t_env *e);

#endif
