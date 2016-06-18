/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_atoi.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/04 11:47:13 by igomez            #+#    #+#             */
/*   Updated: 2014/11/06 15:39:28 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

static int	ft_atoi_rec(int n, const char *str)
{
	if (*str < '0' || *str > '9')
		return (n);
	return (ft_atoi_rec(10 * n + (*str - '0'), str + 1));
}

int			ft_atoi(const char *str)
{
	if (!str)
		return (0);
	while (*str == '\t' || *str == '\n' || *str == '\v' || *str == '\f' ||
			*str == '\r' || *str == ' ')
		str++;
	if (*str == '-')
		return (-ft_atoi_rec(0, str + 1));
	if (*str == '+')
		return (ft_atoi_rec(0, str + 1));
	return (ft_atoi_rec(0, str));
}
