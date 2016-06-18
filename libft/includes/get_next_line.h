/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   get_next_line.h                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/10 10:41:23 by igomez            #+#    #+#             */
/*   Updated: 2015/03/07 15:44:42 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef GET_NEXT_LINE_H
# define GET_NEXT_LINE_H

# define BUFF_SIZE 32

/*
** Get_next_line
*/

typedef struct		s_mem
{
	int				fd_mem;
	char			*data;
}					t_mem;

int					get_next_line(int const fd, char **line);

#endif
