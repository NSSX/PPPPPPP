/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_strmsg.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/26 15:46:52 by igomez            #+#    #+#             */
/*   Updated: 2015/01/28 12:50:11 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

char		*ft_strmsg(int id, t_message *msg_list)
{
	size_t	i;
	char	*msg;

	i = 0;
	if (!msg_list)
		return (NULL);
	while (msg_list[i].msg)
	{
		if (id == msg_list[i++].id)
		{
			msg = msg_list[--i].msg;
			break ;
		}
	}
	if (!msg_list[i].msg)
		msg = "Unknown";
	return (msg);
}
