/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   rtv1_matrix.h                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/05 14:16:05 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 13:47:59 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef RTV1_MATRIX_H
# define RTV1_MATRIX_H

# define M_NEW(M, N)		(t_matrix){M, N, {0}}
# define M_POINT(X, Y, Z)	(t_matrix){4, 1, {X, Y, Z, 1}}
# define M_DIR(X, Y, Z)		(t_matrix){4, 1, {X, Y, Z, 0}}
# define M_IJ(M, I, J)		((M)->mat[I * (M)->n + J])
# define M_PRINT(M)			ft_printmat(M)

/*
** Matrix type, used to represent space tranform (4x4) or vectors (4x1)
*/

typedef struct	s_matrix
{
	int			m;
	int			n;
	double		mat[16];
}				t_matrix;

/*
** Matrix: Calculations
*/

int				mat_prod(t_matrix *buf, const t_matrix *m1, const t_matrix *m2);
double			mat_dot_product(t_matrix *m1, t_matrix *m2);
double			mat_norm(t_matrix *m1);
double			mat_dist(t_matrix *m1, t_matrix *m2);

/*
** Matrix: Space transforms
*/

int				mat_translate(t_matrix *buf, t_matrix *dir, int sign);
int				mat_scale(t_matrix *buf, t_matrix *dir, int sign);
int				mat_normalize(t_matrix *buf);
int				mat_setup(t_matrix *buf, t_matrix *pos, t_matrix *ang,
					t_matrix *scale);
int				mat_cleanup(t_matrix *buf, t_matrix *pos, t_matrix *ang);

/*
** Matrix: Rotations
*/

int				mat_rotate_x(t_matrix *buf, double ang, int sign);
int				mat_rotate_y(t_matrix *buf, double ang, int sign);
int				mat_rotate_z(t_matrix *buf, double ang, int sign);

#endif
